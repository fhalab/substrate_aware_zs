"""
A script for chai structures and also extract the score.
"""

from __future__ import annotations

import os
from glob import glob
from copy import deepcopy

from docko.chai import run_chai

from REVIVAL.preprocess import ZSData
from REVIVAL.util import checkNgen_folder


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
            run_chai(
                label=var,
                seq=seq,
                smiles=self.lib_info["substrate-smiles"]
                + "."
                + ".".join(self.lib_info["cofactor-smiles"]),
                output_dir=checkNgen_folder(
                    os.path.join(self._chai_struct_subdir, var)
                ),
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