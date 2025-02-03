"""
A script for preprocessing the data for the triad and extract energies
"""

import os

from copy import deepcopy
from glob import glob
from tqdm import tqdm

import pandas as pd
import numpy as np

from REVIVAL.global_param import LIB_INFO_DICT
from REVIVAL.preprocess import ZSData
from REVIVAL.util import checkNgen_folder, run_sh_command, get_chain_structure, remove_hetatm


TRIAD_DIR = "/home/bwittmann/triad/triad-2.1.3"


class TriadData(ZSData):
    def __init__(
        self,
        input_csv: str,
        in_structure_dir: str = "data/structure", # data/structure for frompdb or data/structure/docked for docked
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        fit_col_name: str = "fitness",
        triad_dir: str = TRIAD_DIR,
        output_dir: str = "zs/triad",
        withsub: bool = True,
        cleanup: bool = True,
        chain_id: str = "A",
        num_cpus: int = 32,
        rerun: bool = False
    ):
        """
        Run Triad

        For running Triad, there will be subfolders:
            - `struct_prep` for preping pdbs 
            - `mut_file` for .mut files
            - `struct_mut` for generated structure
        For parsing scores, there will be `score` subfolder
        
        Within each, there will be subfolder with info about `withsub` or `cleanup`
            ie `frompdb-cleanup`

        Note that for the docked structure cleanup True or False results in the same structure
        because all HEATM are in chain B LIG and no solvent is present
        """

        super().__init__(
            input_csv=input_csv,
            combo_col_name=combo_col_name,
            var_col_name=var_col_name,
            fit_col_name=fit_col_name,
            withsub=withsub,
        )

        self._in_structure_dir = in_structure_dir
        self._triad_dir = triad_dir
        self._cleanup = cleanup

        if not withsub:
            self._structure_dets = "frompdb"
            if cleanup:
                self._structure_dets += "-cleanup"
        else:
            self._structure_dets = "docked"
            if self._in_structure_dir == "data/structure":
                self._in_structure_dir = os.path.join(self._in_structure_dir, "docked")

        self._output_dir = checkNgen_folder(output_dir)

        self._chain_id = chain_id
        self._num_cpus = num_cpus

        if not rerun and os.path.exists(self.triad_sum_txt):
            print(f"{self.triad_sum_txt} exists and rerun is False")
        
        else:
            self._mut_list = self._generate_mut_file()

            # get the chain structure and run triad preparation
            self._prepare_initial_structure()

            # run triad to generate mutant structures and calculate energies
            self._run_triad()

        # now parse scores
        print(f"Parsing {self.triad_sum_txt} and save to {self.triad_csv}...")

        # extract triad score into dataframe
        # triad_df = self._get_triad_score()

        # # merge to get fit
        # self._triad_df = pd.merge(
        #     self.df[[self._combo_col_name, self._fit_col_name]],
        #     triad_df, 
        #     on =self._combo_col_name,
        #     how="left"
        # )

        # self._triad_df.to_csv(
        #     self.triad_csv, index=False
        # )


    def _generate_mut_file(self) -> None:

        """
        Generate the mut file
        """

        # Loop over variants
        mut_list = []

        for variant in tqdm(self.df[self._var_col_name].tolist()):

            if variant == "WT":
                continue

            # Otherwise, append to mut_list
            mut_list.append(
                "+".join([f"{self._chain_id}_{v[1:]}" for v in variant.split(":")])
                + "\n"
            )

        # check before saving
        # if parent is in the dataframe
        wt_count = self.df[self._var_col_name].eq("WT").sum()

        assert len(mut_list) == len(self.df) - wt_count

        print(f"Generated {self.mut_path}...")

        # Save the mutants
        with open(self.mut_path, "w") as f:
            f.writelines(mut_list)

        return mut_list


    def _prepare_initial_structure(self):
        """Run the proteinProcess.py script to prepare the initial structure."""

        print(f"Getting chain {self._chain_id} structure for {self.instruct_file}...")

        # first get the chain structure
        get_chain_structure(
            input_file_path=self.instruct_file,
            output_file_path=self.triad_instruct_file,
            chain_id=self._chain_id,
        )

        temp_pdb = self.triad_instruct_file.replace(".pdb", "_temp.pdb")

        # clean up hetatm
        if self._cleanup:
            remove_hetatm(input_pdb=self.triad_instruct_file, output_pdb=temp_pdb)
        
            # rename temp to overwrite the self.triad_instruct_file
            os.rename(temp_pdb, self.triad_instruct_file)
        
        # then run traid preparation
        run_sh_command(
            f"{self._triad_dir}/triad.sh {self._triad_dir}/apps/preparation/proteinProcess.py "
            f"-struct {os.path.abspath(self.triad_instruct_file)} "
            f"-crosetta",
            working_dir=os.path.abspath(self.triad_prepstruct_dir)
        )

        # ${TRIAD_DIR}/triad.sh ${TRIAD_DIR}/apps/preparation/proteinProcess.py -struct ${ORIG_PDB_DIR}/${NAME}.pdb -crosetta 


    def _run_triad(self):
        """Run the triad script with the given parameters."""

        run_sh_command(
            f"{self._triad_dir}/tools/openmpi/bin/mpirun -np {str(self._num_cpus)} "
            f"{self._triad_dir}/triad.sh {self._triad_dir}/apps/cleanSequences.py "
            f"-struct {os.path.abspath(self.triad_prepped_file)} "
            f"-rosetta "
            f"-inputSequenceFormat pid "
            f"-inputSequences {os.path.abspath(self.mut_path)} "
            f"-floatNearbyResidues "
            f"-numPDBs={str(self.input_df_length-1)} "
            f"-soft 2>&1 | tee {os.path.abspath(self.triad_sum_txt)}",
            working_dir=os.path.abspath(self.triad_structmut_dir)
            )


    def _get_triad_score(self) -> float:

        """
        A function to load the output of a triad analysis and get a score

        Args:
        - triad_output_file: str, the path to the triad output file
        - WT_combo: str, the wildtype combo
        - num_seqs: int, the number of sequences to load
        """

        # Load the output file
        with open(self.triad_sum_txt) as f:

            # Set some flags for starting analysis
            solutions_started = False
            record_start = False

            # Begin looping over the file
            summary_lines = []
            for line in f:

                # Start looking at data once we hit "solution"
                # if "Solution" in line:
                if "All sequences:" in line:
                    solutions_started = True

                # Once we have "Index" we can start recording the rest
                if solutions_started and "Index" in line:
                    record_start = True

                # Record appropriate data
                if record_start:

                    # Strip the newline and split on whitespace
                    summary_line = line.strip().split()

                    if summary_line[0] == "Average":
                        break
                    else:
                        # Otherwise, append the line
                        summary_lines.append(summary_line)

        # Build the dataframe with col ['Index', 'Tags', 'Score', 'Seq', 'Muts']
        all_results = pd.DataFrame(summary_lines[1:], columns=summary_lines[0])
        all_results["Triad_score"] = all_results["Score"].astype(float)

        wt_chars = self.parent_aa
        reconstructed_combos = [
            "".join(
                [char if char != "-" else wt_chars[i] for i, char in enumerate(seq)]
            )
            for seq in all_results.Seq.values
        ]
        all_results[self._combo_col_name] = reconstructed_combos

        # Get the order
        all_results["Triad_rank"] = np.arange(1, len(all_results) + 1)

        return all_results[[self._combo_col_name, "Triad_score", "Triad_rank"]]


    @property
    def triad_mut_dir(self) -> str:
        """
        A property for the triad mut directory
        """
        return checkNgen_folder(
            os.path.join(self._output_dir, "mut_file", self._structure_dets)
        )

    @property
    def triad_prepstruct_dir(self) -> str:
        """
        A property for the triad struct directory
        """
        return checkNgen_folder(
            os.path.join(self._output_dir, "struct_prep", self._structure_dets, self.zs_struct_name)
        )
    
    @property
    def instruct_file(self) -> str:

        """
        PDB file path to the triad pdb file
        """

        struct_ext = "pdb" if "frompdb" in self._structure_dets else "cif"

        return os.path.join(
            self._in_structure_dir,
            f"{self.zs_struct_name}.{struct_ext}",
        )

    @property
    def triad_instruct_file(self) -> str:

        """
        PDB file path to the triad pdb file
        """

        return os.path.join(self.triad_prepstruct_dir, f"{self.zs_struct_name}.pdb")

    @property
    def triad_prepped_file(self) -> str:
        """
        PDB file path to the triad pdb file
        """

        return self.triad_instruct_file.replace(".pdb", "_prepared.pdb")

    @property
    def triad_structmut_dir(self) -> str:
        """
        A property for the triad output directory
        """
        return checkNgen_folder(
            os.path.join(self._output_dir, "struct_mut", self._structure_dets, self.lib_name)
        )
    
    @property
    def triad_sum_txt(self) -> str:
        """
        A property for the triad summary txt file path
        """
        return os.path.join(self.triad_structmut_dir, f"{self.lib_name}.txt")

    @property
    def prefixes(self) -> list:
        """
        A property for the prefixes
        """
        return [
            f"{self._chain_id}_{pos}"
            for pos in LIB_INFO_DICT[self.lib_name]["positions"].values()
        ]

    @property
    def mut_path(self) -> str:
        """
        A property for the mut file path
        """
        return os.path.join(self.triad_mut_dir, f"{self.lib_name}.mut")

    @property
    def mut_list(self) -> list:
        """
        A property for the mutation encodings
        """
        return self._mut_list

    @property
    def triad_score_dir(self) -> str:
        """Folder for csv score files"""
        return checkNgen_folder(
            os.path.join(self._output_dir, "score", self._structure_dets)
        )

    @property
    def triad_csv(self) -> str:
        """
        A property for the triad csv
        """
        return os.path.join(self.triad_score_dir, f"{self.lib_name}.csv")

    @property
    def triad_df(self) -> pd.DataFrame:
        """
        A property for the triad dataframe
        """
        return self._triad_df


def run_traid(
    pattern = "data/meta/not_scaled/*.csv",
    in_structure_dir = "data/structure",
    kwargs: dict = {}
):
    """
    Run the triad gen mut file function for all libraries

    Args:
    - pattern: str | list: the pattern for the input csv files
    """

    if isinstance(pattern, str):
        lib_list = sorted(glob(pattern))
    else:
        lib_list = deepcopy(pattern)

    for lib in tqdm(lib_list):
        print(f"Running triad for {lib}...")
        TriadData(input_csv=lib, in_structure_dir=in_structure_dir, **kwargs)
