"""
A script for preprocessing the data for the triad and extract energies
"""

import re
import os

from copy import deepcopy
from glob import glob

import pandas as pd
import numpy as np
from itertools import product

from REVIVAL.global_param import LIB_INFO_DICT
from REVIVAL.preprocess import ZSData
from REVIVAL.util import checkNgen_folder, get_file_name


class TriadData(ZSData):
    def __init__(self,
        input_csv: str,
        scale_fit: str = "parent",
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        mut_col_name: str = "mut",
        pos_col_name: str = "pos",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        seq_dir: str = "data/seq",
        zs_dir: str = "zs",
        triad_dir: str = "triad",
        triad_mut_dir: str = "mut_file",
        chain_id: str = "A",
    ):

        super().__init__(
            input_csv,
            scale_fit,
            combo_col_name,
            var_col_name,
            mut_col_name,
            pos_col_name,
            seq_col_name,
            fit_col_name,
            seq_dir,
            zs_dir
        )

        self._triad_dir = checkNgen_folder(os.path.join(self._zs_dir, triad_dir))
        self._triad_mut_dir = checkNgen_folder(os.path.join(self._triad_dir, triad_mut_dir))
        self._chain_id = chain_id

        self._mut_list = self._generate_mut_file()

    def _generate_mut_file(self) -> None:

        """
        Generate the mut file
        """

        # Loop over variants
        mut_list = []

        for variant in self.df[self._combo_col_name].tolist():

            # Loop over each character in the variant
            mut_char_list = []
            for j, (var_char, wt_char) in enumerate(zip(variant, self.parent_aa)):

                # If the var_char does not equal the wt_char, append
                if var_char != wt_char:
                    mut_char_list.append(self.prefixes[j] + var_char)

            # If the mut_char_list has no entries, continue (this is wild type)
            if len(mut_char_list) == 0:
                continue

            # Otherwise, append to mut_list
            else:
                mut_list.append("+".join(mut_char_list) + "\n")

        # check before saving
        assert len(mut_list) == self.df_length - 1

        print(f"Generated {self.mut_path}...")

        # Save the mutants
        with open(self.mut_path, "w") as f:
            f.writelines(mut_list)

        return mut_list

    @property
    def prefixes(self) -> list:
        """
        A property for the prefixes
        """
        return [
            f"{self._chain_id}_{pos}" for pos in LIB_INFO_DICT[self.lib_name]["positions"].values()
        ]

    @property
    def mut_path(self) -> str:
        """
        A property for the mut file path
        """
        return os.path.join(self._triad_mut_dir, f"{self.lib_name}.mut")

    @property
    def mut_list(self) -> list:
        """
        A property for the mutation encodings
        """
        return self._mut_list


class TriadResults(ZSData):
    """
    A class for parsing the triad results
    """

    def __init__(
        self,
        input_csv: str,
        scale_fit: str = "parent",
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        mut_col_name: str = "mut",
        pos_col_name: str = "pos",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        seq_dir: str = "data/seq",
        zs_dir: str = "zs",
        triad_dir: str = "triad",
        triad_rawouput_dir: str = "raw_output",
        triad_processed_dir: str = "processed_output",
        chain_id: str = "A",
    ) -> None:

        """
        Args:
        """

        super().__init__(
                        input_csv,
            scale_fit,
            combo_col_name,
            var_col_name,
            mut_col_name,
            pos_col_name,
            seq_col_name,
            fit_col_name,
            seq_dir,
            zs_dir
        )

        self._triad_rawouput_dir = checkNgen_folder(os.path.join(zs_dir, triad_dir, triad_rawouput_dir))
        self._triad_processed_dir = checkNgen_folder(os.path.join(zs_dir, triad_dir, triad_processed_dir))

        print(f"Parsing {self.triad_txt} and save to {self.triad_csv}...")

        # extract triad score into dataframe
        self._triad_df = self._get_triad_score()


    def _get_triad_score(self) -> float:

        """
        A function to load the output of a triad analysis and get a score

        Args:
        - triad_output_file: str, the path to the triad output file
        - WT_combo: str, the wildtype combo
        - num_seqs: int, the number of sequences to load
        """

        # Load the output file
        with open(self.triad_txt) as f:

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
        all_results["AAs"] = reconstructed_combos

        # Get the order
        all_results["Triad_rank"] = np.arange(1, len(all_results) + 1)

        return all_results[["AAs", "Triad_score", "Triad_rank"]].to_csv(self.triad_csv, index=False)

    @property
    def triad_txt(self) -> str:
        """
        A property for the triad txt
        """
        return os.path.join(self._triad_rawouput_dir, f"{self.lib_name}.txt")
    
    @property
    def triad_csv(self) -> str:
        """
        A property for the triad csv
        """
        return os.path.join(self._triad_processed_dir, f"{self.lib_name}.csv")

    @property
    def triad_df(self) -> pd.DataFrame:
        """
        A property for the triad dataframe
        """
        return self._triad_df


def run_traid_gen_mut_file(
    pattern: str | list = "data/meta/scale2parent/*.csv",
    kwargs: dict = {}
    ):
    """
    Run the triad gen mut file function for all libraries
    
    Args:
    - pattern: str | list: the pattern for the input csv files
    """

    if isinstance(pattern, str):
        lib_list = glob(pattern)
    else:
        lib_list = deepcopy(pattern)

    for lib in lib_list:
        print(f"Running triad gen mut file for {lib}...")
        TriadData(input_csv=lib, **kwargs)


def run_parse_triad_results(
    pattern: str | list = "data/meta/scale2parent/*.csv",
    kwargs: dict = {}
    ):

    """
    Run the parse triad results function for all libraries
    
    Args:
    """
    if isinstance(pattern, str):
        lib_list = glob(pattern)
    else:
        lib_list = deepcopy(pattern)

    for lib in lib_list:
        print(f"Running parse triad results for {lib}...")
        TriadResults(input_csv=lib, **kwargs)
