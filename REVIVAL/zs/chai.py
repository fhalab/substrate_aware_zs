"""
A script for chai structures and also extract the score.
"""

from __future__ import annotations

import os
from glob import glob
from copy import deepcopy

import numpy as np
import pandas as pd

from docko.chai import run_chai

from REVIVAL.preprocess import ZSData
from REVIVAL.util import checkNgen_folder, get_file_name


class ChaiData(ZSData):
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
        chai_dir: str = "chai",
        chai_struct_dir: str = "mut_structure",
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
            zs_dir,
        )

        self._chai_dir = checkNgen_folder(os.path.join(self._zs_dir, chai_dir))
        self._chai_struct_dir = checkNgen_folder(
            os.path.join(self._chai_dir, chai_struct_dir)
        )
        self._chai_struct_subdir = checkNgen_folder(
            os.path.join(self._chai_struct_dir, self.lib_name)
        )

        self._gen_chai_structure()

    def _gen_chai_structure(self):
        """
        A method to generate the chai structure for each variant.
        """

        for (
            var,
            seq,
        ) in self.df[[self._var_col_name, self._seq_col_name]].values:

            print(f"Running chai for {var} saving to {self._chai_struct_subdir}...")

            run_chai(
                label=var,
                seq=seq,
                smiles=self.lib_info["substrate-smiles"]
                + "."
                + ".".join(self.lib_info["cofactor-smiles"]),
                output_dir=self._chai_struct_subdir,
            )


def run_gen_chai_structure(
    pattern: str | list = "data/meta/scale2parent/*.csv", kwargs: dict = {}
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
        ChaiData(input_csv=lib, **kwargs)


def parse_chai_scores(mut_structure_dir: str, output_dir: str = "zs/chai/output"):

    """
    A function for going through the subfolder and getting the chai scores
    to generate a dataframe with the following columns:
        - var: The mutation, ie I165A:I183A:Y301V
        - rep: The replicate number
        - aggregate_score
        - ptm
        - iptm
        - chain_ptm_A
        - chain_ptm_B
        - chain_iptm_AA
        - chain_iptm_AB
        - chain_iptm_BA
        - chain_iptm_BB
        - has_inter_chain_clashes

    Args:
    - input_dir, str: The path to the folder containing the chai score
        ie zs/chai/mut_structure/PfTrpB-4bromo
    - output_dir, str: The path to the folder to save the dataframe to
        ie zs/chai/output
    """

    output_dir = checkNgen_folder(output_dir)
    lib_name = os.path.basename(mut_structure_dir)

    # init dataframe
    df = pd.DataFrame(
        columns=[
            "var",
            "rep",
            "aggregate_score",
            "ptm",
            "iptm",
            "chain_ptm_A",
            "chain_ptm_B",
            "chain_iptm_AB",
            "chain_iptm_BA",
            "has_inter_chain_clashes",
        ]
    )

    for subfolder in glob(f"{mut_structure_dir}/*"):
        var = os.path.basename(subfolder)
        
        for rep_npz in glob(f"{subfolder}/*.npz"):

            npz = np.load(rep_npz)

            df = df._append(
                {
                    "var": var,
                    "rep": get_file_name(rep_npz).split("_")[-1],
                    "aggregate_score": npz["aggregate_score"][0],
                    "ptm": npz["ptm"][0],
                    "iptm": npz["iptm"][0],
                    "chain_ptm_A": npz["per_chain_ptm"][0][0],
                    "chain_ptm_B": npz["per_chain_ptm"][0][1],
                    "chain_iptm_AB": npz["per_chain_pair_iptm"][0][0, 1],
                    "chain_iptm_BA": npz["per_chain_pair_iptm"][0][1, 0],
                    "has_inter_chain_clashes": npz["has_inter_chain_clashes"][0],
                },
                ignore_index=True,
            )

    df.to_csv(f"{output_dir}/{lib_name}.csv", index=False)
    print(f"Saved chai scores for {lib_name} to {output_dir}/{lib_name}.csv")


def parse_all_chai_scores(chai_struct_dir: str = "zs/chai/mut_structure"):
    """
    A function to parse all the chai scores for all libraries
    """
    
    for lib in glob(f"{chai_struct_dir}/*"):
        print(f"Parsing chai scores for {lib}...")
        parse_chai_scores(lib)