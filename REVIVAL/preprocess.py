"""
Preprocess the input data
"""

from __future__ import annotations

import os
from ast import literal_eval
from glob import glob
from tqdm import tqdm
from copy import deepcopy

import numpy as np
import pandas as pd

from REVIVAL.util import checkNgen_folder, get_file_name, read_parent_fasta
from REVIVAL.global_param import LIB_INFO_DICT


class LibData:
    """
    A parent class to get the library information
    """

    def __init__(
        self,
        input_csv: str,
        scale_fit: str,
        combo_col_name: str = "AAs",
        mut_col_name: str = "muts",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        seq_dir: str = "data/seq",
    ) -> None:
        """
        Args:
        - input_csv, str: path to the input csv file,
            ie. data/lib/PfTrpB-4bromo.csv for preprocessed
                data/meta/scale2max/PfTrpB-4bromo.csv for scaled to max = 1
        - scale_fit, str: ways to scale the fitness
            'parent' means the parent fitness = 1
            'max' means max fitness = 1
        - combo_col_name, str: the column name for the mutated amino acids
        - mut_col_name, str: the column name for the mutations
        - seq_col_name, str: the column name for the full sequence
        - fit_col_name, str: the column name for the fitness
        - seq_dir, str: the directory for the parent sequence fasta files
        """

        self._input_csv = input_csv
        self._scale_fit = scale_fit

        self._combo_col_name = combo_col_name
        self._mut_col_name = mut_col_name
        self._seq_col_name = seq_col_name
        self._fit_col_name = fit_col_name
        self._seq_dir = seq_dir

    @property
    def lib_name(self) -> dict:
        """Return the library name"""
        return get_file_name(self._input_csv)

    @property
    def lib_info(self) -> dict:
        """Return the library information"""
        return LIB_INFO_DICT[self.lib_name]

    @property
    def protein_name(self) -> str:
        """
        Returns the protein name
        """
        return self.lib_info["enzyme"]

    @property
    def parent_aa(self) -> str:
        """Return the parent amino acid"""
        return "".join(list(self.lib_info[self._combo_col_name].values()))

    @property
    def n_site(self) -> int:
        """Return the number of sites"""
        return len(self.lib_info["positions"])

    @property
    def scale_type(self) -> str:
        """Return the scale type"""
        if self._scale_fit in ["max", "parent"]:
            return f"scale2{self._scale_fit}"
        else:
            return "not_scaled"

    @property
    def input_df(self) -> pd.DataFrame:
        """Return the input dataframe"""
        return pd.read_csv(self._input_csv)

    @property
    def df_length(self):
        return len(self.input_df)

    @property
    def split_aa_cols(self) -> list:
        """Return the columns for the split amino acids"""
        return [f"AA{str(i)}" for i in self.lib_info["positions"].keys()]

    @property
    def seq_file(self) -> str:
        """Return the path to the parent sequence fasta file"""
        return os.path.join(self._seq_dir, f"{self.protein_name}.fasta")

    @property
    def parent_seq(self) -> str:
        """Return the parent sequence"""
        return read_parent_fasta(self.seq_file)


class ProcessData(LibData):
    """
    A parent class to process the data
    """

    def __init__(
        self,
        input_csv: str,
        scale_fit: str = "parent",
        combo_col_name: str = "AAs",
        mut_col_name: str = "muts",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        seq_dir: str = "data/seq",
        output_dir: str = "data/meta",
    ) -> None:

        """
        Args:
        - input_csv, str: path to the input csv file
        - scale_fit, str: ways to scale the fitness
            'parent' means the parent fitness = 1
            'max' means max fitness = 1
        """

        super().__init__(
            input_csv,
            scale_fit,
            combo_col_name,
            mut_col_name,
            seq_col_name,
            fit_col_name,
            seq_dir,
        )

        self._output_subdir = checkNgen_folder(
            os.path.join(output_dir, self.scale_type)
        )

        print(f"Processing {self._input_csv} with {self._scale_fit}...")

        # save the output csv
        print(f"Saving to {self.output_csv} ...")

    def _split_aa(self, df: pd.DataFrame) -> pd.DataFrame:

        """
        Split the amino acids into individual columns

        Args:
        - df, pd.DataFrame: the input dataframe

        Returns:
        - pd.DataFrame: the split dataframe
        """

        df_split = df.copy()

        df_split[self.split_aa_cols] = df_split[self._combo_col_name].apply(
            lambda x: pd.Series(list(x))
        )

        return df_split.copy()

    def _convert_muts(self, muts: str) -> str:

        """
        Convert the variants sequence
        to the form of parentaa1loc1mutaa1:parentaa2loc2mutaa2
        for mut number count

        Args:
        - muts, str: the variants sequence

        Returns:
        - str: the converted sequence
        """

        mut_seq = ""
        mut_num = 0

        for i, (mut, wt) in enumerate(zip(muts, self.parent_aa)):
            # get wt aa + mut position from ie
            # "TrpB4": {
            #     "positions": {1: 183, 2: 184, 3: 227, 4: 228},
            #     "codons": {1: "GTG", 2: "TTC", 3: "GTG", 4: "AGC"},
            #     self._combo_col_name: {1: "V", 2: "F", 3: "V", 4: "S"}
            # }

            wt_aa_pos = LIB_INFO_DICT[self.lib_name][self._combo_col_name][i + 1] + str(
                LIB_INFO_DICT[self.lib_name]["positions"][i + 1]
            )

            if mut != wt:
                mut_num += 1
                if mut_num != 1:
                    mut_seq += ":" + wt_aa_pos + mut
                else:
                    mut_seq += wt_aa_pos + mut
        return mut_seq

    def _append_mut(self, df: pd.DataFrame) -> pd.DataFrame:

        """
        Apply the convert_muts function to the dataframe

        Args:
        - df, pd.DataFrame: the input dataframe

        Returns:
        - pd.DataFrame: the appended dataframe
        """

        df_appended = df.copy()
        df_appended.loc[:, self._mut_col_name] = df_appended.apply(
            lambda x: self._convert_muts(x[self._combo_col_name]),
            axis=1,
        )

        df_appended[self._mut_col_name] = df_appended[self._mut_col_name].replace(
            "", "WT"
        )

        # add mut number
        df_appended["n_mut"] = df_appended[self._mut_col_name].str.split(":").str.len()

        # change WT n_mut to 0
        df_appended.loc[df_appended[self._mut_col_name] == "WT", "n_mut"] = 0

        return df_appended.copy()

    def _mut2seq(self, muts: str) -> str:

        """
        Apply the mutations to the parent sequence

        Args:
        - muts, str: the mutations that are 1 indexed
            ie. C1A:A3D

        Returns:
        - str: the mutated sequence
        """

        seq_list = list(self.parent_seq)

        if muts == "WT":
            return self.parent_seq

        else:
            for mut in muts.split(":"):

                # first character is the parent amino acid and last is the mutated amino acid
                wt_aa, mut_aa = mut[0], mut[-1]

                # Extract position (1-indexed, convert to 0-indexed)
                mut_pos = int(mut[1:-1]) - 1

                # Assert that the parent sequence has the expected amino acid at the position
                assert (
                    seq_list[mut_pos] == wt_aa
                ), f"Mismatch at position {mut_pos + 1}: expected {wt_aa}, found {seq_list[mut_pos]}"

                # Replace the parent amino acid with the mutated amino acid at the specified position
                seq_list[mut_pos] = mut_aa

            return "".join(seq_list)

    def _process(self) -> pd.DataFrame:

        """
        Process the input csv
        """

        df = self.input_df.copy()

        # append muts column for none SSM data
        if self._mut_col_name not in df.columns and self._combo_col_name in df.columns:
            df = self._append_mut(df).copy()

        # split the amino acids for SSM data
        if "AA1" not in df.columns and self._combo_col_name in df.columns:
            df = self._split_aa(df).copy()

        # add full seq from fasta file by modifying self.parent_seq with the mutations
        if self._seq_col_name not in df.columns and self._combo_col_name in df.columns:
            df.loc[:, self._seq_col_name] = df[self._mut_col_name].apply(lambda x: self._mut2seq(x))

        # add col for enzyme name, substrate, cofactor, and their smile strings if relevant
        for col in ["enzyme", "substrate", "substrate-smiles", "cofactor", "cofactor-smiles"]:
            if col in self.lib_info:
                df[col] = self.lib_info[col]

        # save the output csv
        df.to_csv(self.output_csv, index=False)

    @property
    def output_csv(self) -> str:

        """Return the path to the output csv"""

        return os.path.join(self._output_subdir, os.path.basename(self._input_csv))


def preprocess_all(
    input_pattern: str = "data/lib/*.csv",
    scale_fit: str = "parent",
    combo_col_name: str = "AAs",
    mut_col_name: str = "muts",
    seq_col_name: str = "seq",
    fit_col_name: str = "fitness",
    seq_dir: str = "data/seq",
    output_dir: str = "data/meta",
) -> None:

    """
    Preprocess all the csv files in the input directory

    Args:
    - input_pattern, str: the input directory
    - scale_fit, str: ways to scale the fitness
        'parent' means the parent fitness = 1
        'max' means max fitness = 1
    - combo_col_name, str: the column name for the mutated amino acids
    - mut_col_name, str: the column name for the mutations
    - seq_col_name, str: the column name for the full sequence
    - fit_col_name, str: the column name for the fitness
    - seq_dir, str: the directory for the parent sequence fasta files
    - output_dir, str: the output directory

    """

    for input_csv in tqdm(glob(input_pattern)):

        ProcessData(
            input_csv=input_csv,
            scale_fit=scale_fit,
            combo_col_name=combo_col_name,
            mut_col_name=mut_col_name,
            seq_col_name=seq_col_name,
            fit_col_name=fit_col_name,
            seq_dir=seq_dir,
            output_dir=output_dir,
        )._process()