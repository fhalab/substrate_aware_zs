"""
A script for combining the results of multiple zs.
"""

import os
import re
import sys
from glob import glob
from tqdm import tqdm
from copy import deepcopy

import pandas as pd
import numpy as np


from REVIVAL.preprocess import ZSData
from REVIVAL.util import checkNgen_folder


# Regex pattern to ignore columns ending with "_0" or "_1"
IGNORE_PATTERN = re.compile(r".*_(\d+)$")


class ZSComb(ZSData):
    """
    Combine multiple zs scores
    """

    def __init__(
        self,
        input_csv: str,
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        fit_col_name: str = "fitness",
        zs_dir: str = "zs",
        comb_dir: str = "comb",
        zs_subdir_list: list = [
            "ev/*.csv",
            "esm/output/*.csv",
            "flowsite/scores/pocket_def_residues_model2-substrate-cofactor/*.csv",
            "ligandmpnn/scores/autoregressive_20/*.csv",
            "vol/*.csv",  #
            "esmif/*/score/*.csv",
            "coves/*/output/100_processed/*.csv",
            "triad/score/*/*.csv",
            "af3/score_*/*.csv",  # need to look at apo
            "chai/score_*/*.csv",  # need to look at apo
            "bonddist/*/*/*.csv",
            "hydro/**/*.csv",
            "plip/*/*/*.csv",
            "vina/*/*/*/*.csv",  # add af3 chai sub dock from smiles
        ],
    ):

        super().__init__(
            input_csv=input_csv,
            combo_col_name=combo_col_name,
            var_col_name=var_col_name,
            fit_col_name=fit_col_name,
            zs_dir=zs_dir,
        )

        self._comb_dir = checkNgen_folder(os.path.join(zs_dir, comb_dir))
        self._zs_subdir_list = deepcopy(zs_subdir_list)

        self._comb_zs = self._comb_zs()
        self._comb_zs.to_csv(self.zs_comb_path, index=False)

    def _set_pattern_property(self):
        """
        For each item in self.zs_opts, set the pattern property
        using the pattern in self._zs_subdir_list

        ie self._ev_pattern = "ev/*.csv"
        """
        for opt in self.zs_opts:
            setattr(
                self, f"_{opt}_pattern", self._zs_subdir_list[self.zs_opts.index(opt)]
            )

    def _append_complex_zs_csv(self, pattern):
        """
        Append all type of unique option for a complex zs
        """

        # Store all DataFrames for merging
        combined_df = pd.DataFrame()

        # Use glob to find matching CSV files
        csv_files = glob(pattern, recursive=True)
        

        # Extract the unique portions of each path
        for csv_path in csv_files:

            split_csv_path = csv_path.split("/")

            # Replace `/` with `-` to create a unique identifier
            unique_name = "-".join(split_csv_path[2:-1])

            # to prevent duplicate names from af3 and chai outputs
            if split_csv_path[1] in ["af3", "chai"]:
                unique_name += f"_{split_csv_path[1]}"

            print(f"Combining {csv_path} with {unique_name}...")

            # Load CSV file
            df = pd.read_csv(csv_path)

            # Remove columns that match IGNORE_PATTERN
            df = df.loc[:, ~df.columns.str.match(IGNORE_PATTERN)]

            # Rename columns by appending unique_combo (unless in EXCLUDE_COLUMNS)
            df = df.rename(
                columns={
                    col: f"{col}_{unique_name}" if col not in self.common_cols else col
                    for col in df.columns
                }
            )

            combine_col = [
                c
                for c in self.common_cols
                if c in df.columns and c in combined_df.columns
            ]

            # Store processed DataFrame
            combined_df = (
                pd.merge(combined_df, df, on=combine_col, how="outer")
                if not combined_df.empty
                else df
            )

        return combined_df

    def _comb_zs(self) -> pd.DataFrame:

        """
        Combine the zs scores
        """

        df = self.df[self.common_cols].copy()

        # add hamming distance first
        df["hd"] = -1 * df["n_mut"]

        # first get the simple ones
        for zs_path in [
            os.path.join(self._zs_dir, f)
            for f in self._zs_subdir_list
            if f.count("*") == 1
        ]:

            zs_path = zs_path.replace("*.csv", f"{self.lib_name}.csv")

            print(f"Combining {zs_path}...")

            zs_df = pd.read_csv(zs_path)

            simple_common_cols = list(set(df.columns) & set(zs_df.columns))
            print(f"on {simple_common_cols}...")
            df = pd.merge(df, zs_df, on=simple_common_cols, how="outer")

        # then get the complex ones
        for zs_path in [
            os.path.join(self._zs_dir, f)
            for f in self._zs_subdir_list
            if f.count("*") > 1
        ]:

            zs_path = zs_path.replace("*.csv", f"{self.lib_name}.csv")

            print(f"Combining {zs_path}...")
            complex_zs_df = self._append_complex_zs_csv(zs_path)
            complex_common_cols = list(set(df.columns) & set(complex_zs_df.columns))
            print(f"on {complex_common_cols}...")
            df = pd.merge(df, complex_zs_df, on=complex_common_cols, how="outer")

        return df.copy()

    @property
    def zs_opts(self) -> list:
        """
        Return the list of zs subdirectories
        """
        return deepcopy([f.split("/")[0] for f in self._zs_subdir_list])

    @property
    def zs_comb_path(self) -> str:
        """
        Return the path to the combined zs
        """
        return os.path.join(self._comb_dir, self.lib_name + ".csv")

    @property
    def complex_zs_opts(self) -> list:
        """
        Find those has more than one * in the path
        """
        return [
            os.path.join(self._zs_dir, f)
            for f in self._zs_subdir_list
            if f.count("*") > 1
        ]

    @property
    def common_cols(self) -> list:
        """
        Return the list of common columns
        """
        common_cols = [self._combo_col_name, self._var_col_name, self._fit_col_name]
        add_cols = ["selectivity", "n_mut", "enzyme"]
        for c in add_cols:
            if c in self.df.columns:
                common_cols.append(c)
        return common_cols


def run_all_combzs(
    pattern: str = "data/meta/not_scaled/*",
    combo_col_name: str = "AAs",
    var_col_name: str = "var",
    fit_col_name: str = "fitness",
    zs_dir: str = "zs",
    comb_dir: str = "comb",
):

    """
    Combine all scores for all datasets
    """
    if isinstance(pattern, str):
        path_list = sorted(glob(pattern))
    else:
        path_list = deepcopy(pattern)

    for p in tqdm(path_list):

        print(f"Running zs comb for {p}...")

        ZSComb(
            input_csv=p,
            combo_col_name=combo_col_name,
            var_col_name=var_col_name,
            fit_col_name=fit_col_name,
            zs_dir=zs_dir,
            comb_dir=comb_dir,
        )