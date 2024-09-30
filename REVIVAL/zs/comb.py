"""
A script for combining the results of multiple zs.
"""

import os
import sys
from glob import glob
from copy import deepcopy

import pandas as pd
import numpy as np


from REVIVAL.preprocess import ZSData
from REVIVAL.util import checkNgen_folder

class ZSComb(ZSData):
    """
    Combine multiple zs scores
    """

    def __init__(
        self,
        input_csv: str,
        scale_fit: str,
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        mut_col_name: str = "mut",
        pos_col_name: str = "pos",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        seq_dir: str = "data/seq",
        zs_dir: str = "zs",
        comb_dir: str = "comb",
        zs_subdir_list: list = ["esm/output", "esmif/output"],
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

        self._comb_dir = checkNgen_folder(os.path.join(zs_dir, comb_dir))
        self._zs_subdir_list = zs_subdir_list

        self._comb_zs()

    def _comb_zs(self) -> pd.DataFrame:

        """
        Combine the zs scores
        """

        df = self.input_df.copy()
        for zs_path in self.zs_paths:
            zs_df = pd.read_csv(zs_path)
            common_cols = list(set(df.columns) & set(zs_df.columns))
            df = pd.merge(
                df,
                zs_df,
                on=common_cols,
            )

        df.to_csv(self.zs_comb_path, index=False)

        return df.copy()
        

    @property
    def zs_opts(self) -> list:
        """
        Return the list of zs subdirectories
        """
        return deepcopy([f.split("/")[0] for f in self._zs_subdir_list])

    @property
    def zs_paths(self) -> list:
        """
        Return the list of zs paths
        """
        return deepcopy([os.path.join(self._zs_dir, f, self.lib_name + ".csv") for f in self._zs_subdir_list])

    @property
    def zs_comb_path(self) -> str:
        """
        Return the path to the combined zs
        """
        return os.path.join(self._comb_dir, self.lib_name + ".csv")


def run_all_combzs(
    pattern: str = "data/meta/scale2parent/*",
    scale_fit: str = "parent",
    combo_col_name: str = "AAs",
    var_col_name: str = "var",
    mut_col_name: str = "mut",
    pos_col_name: str = "pos",
    seq_col_name: str = "seq",
    fit_col_name: str = "fitness",
    seq_dir: str = "data/seq",
    zs_dir: str = "zs",
    comb_dir: str = "comb",
    zs_subdir_list: list = ["esm/output"],
):

    """
    Combine all scores for all datasets
    """
    if isinstance(pattern, str):
        path_list = glob(pattern)
    else:
        path_list = deepcopy(pattern)

    for p in path_list:

        print(f"Running ESM for {p}...")

        ZSComb(
            input_csv=p,
            scale_fit=scale_fit,
            combo_col_name=combo_col_name,
            var_col_name=var_col_name,
            mut_col_name=mut_col_name,
            pos_col_name=pos_col_name,
            seq_col_name=seq_col_name,
            fit_col_name=fit_col_name,
            seq_dir=seq_dir,
            zs_dir=zs_dir,
            comb_dir=comb_dir,
            zs_subdir_list=zs_subdir_list
        )