"""
A script for generating fasta files for ESMIF.
"""

import os
from glob import glob
from copy import deepcopy

from REVIVAL.preprocess import ZSData
from REVIVAL.util import convert_cif_to_pdb, checkNgen_folder


class ESMIFData(ZSData):

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
        structure_dir: str = "data/structure",
        zs_dir: str = "zs",
        esmif_dir: str = "esmif",
        esmif_mut_dir: str = "mut_file",
        esmif_ipdb_dir: str = "input_pdb",
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
            structure_dir,
            zs_dir
        )

        self._esmif_dir = checkNgen_folder(os.path.join(self._zs_dir, esmif_dir))
        self._esmif_mut_dir = checkNgen_folder(os.path.join(self._esmif_dir, esmif_mut_dir))
        self._esmif_ipdb_dir = checkNgen_folder(os.path.join(self._esmif_dir, esmif_ipdb_dir))
        self._chain_id = chain_id

        # process input pdb
        self._process_input_pdb()


    def _process_input_pdb(self) -> None:
        """
        A function for processing input pdb

        Args:
        - input_pdb: str, path to input pdb
        """

        if ".cif" in self.structure_file:
            print(f"Converting {self.structure_file} to {self.esmif_ipdb_file}...")
            convert_cif_to_pdb(self.structure_file, self.esmif_ipdb_file)
        else:
            # copy the input pdb to the esmif input pdb directory
            print(f"Copying {self.structure_file} to {self.esmif_ipdb_file}...")
            os.system(f"cp {self.structure_file} {self.esmif_ipdb_file}")


    def _mut_csv2fasta(self) -> None:
        """
        A function for converting mutation csv to fasta
        """

        # make sure the seq and mut names are in the csv
        for col in [self._var_col_name, self._seq_col_name]:
            if col not in self.df.columns:
                raise ValueError(f"{col} column not found")

        print(f"Writing to {self.esmif_mut_file}...")

        # TODO: dealt with seq not same length
        with open(self.esmif_mut_file, "w") as f:
            for mut, seq in zip(self.df[self._var_col_name].values, self.df[self._seq_col_name].values):
                f.write(f">{mut}\n{seq}\n")
    
    @property
    def esmif_ipdb_file(self) -> str:
        """
        PDB file path to the esmif pdb file
        """
        return os.path.join(self._esmif_ipdb_dir, f"{self.lib_name}.pdb")

    @property
    def esmif_mut_file(self) -> str:
        """
        Mutation file path to the esmif mutation file
        """
        return os.path.join(self._esmif_mut_dir, f"{self.lib_name}.fasta")


def run_esmif_input_file(
    pattern: str | list = "data/meta/scale2parent/*.csv",
    kwargs: dict = {}
    ):
    """
    Run the esmif input file function for all libraries
    
    Args:
    - pattern: str | list: the pattern for the input csv files
    """

    if isinstance(pattern, str):
        lib_list = glob(pattern)
    else:
        lib_list = deepcopy(pattern)

    for lib in lib_list:
        print(f"Running gen mut file for {lib}...")
        ESMIFData(input_csv=lib, **kwargs)._mut_csv2fasta()
