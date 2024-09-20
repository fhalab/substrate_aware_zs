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
        sub_folder = checkNgen_folder(os.path.join(self._triad_dir, "mut_file"))
        return os.path.join(sub_folder, f"{self.lib_name}.mut")

    @property
    def mut_list(self) -> list:
        """
        A property for the mutation encodings
        """
        return self._mut_list


def run_traid_gen_mut_file(
    pattern: str | list = "data/meta/scale2parent/*.csv"
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
        TriadData(input_csv=lib)


