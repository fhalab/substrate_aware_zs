"""
Preprocess the input data
"""

from __future__ import annotations

import os
from glob import glob
from tqdm import tqdm
from copy import deepcopy

import numpy as np
import pandas as pd

from substrate_aware.global_param import LIB_INFO_DICT, APPEND_INFO_COLS
from substrate_aware.chem_helper import canonicalize_smiles
from substrate_aware.util import (
    checkNgen_folder,
    get_file_name,
    read_parent_fasta,
    er2ee,
    get_chain_structure,
    remove_hetatm,
    remove_residues_from_pdb,
)


class LibData:
    """
    A parent class to get the library information
    """

    def __init__(
        self,
        input_csv: str,
        scale_fit: str = "",
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        selectivity_col_name: str = "er",
        protein_name: str = "",
        seq_dir: str = "data/seq",
        structure_dir: str = "data/structure",
        mut_fasta_dir: str = "data/mut_fasta",
    ) -> None:

        """
        Args:
        - input_csv, str: path to the input csv file,
            ie. data/lib/PfTrpB-4bromo.csv for preprocessed
                data/meta/scale2max/PfTrpB-4bromo.csv for scaled to max = 1
        - scale_fit, str: ways to scale the fitness
            'parent' means the parent fitness = 1
            'max' means max fitness = 1
            '' or 'not_scaled' means no scaling, default
        - combo_col_name, str: the column name for the mutated amino acids
            ie 'VDGV'
        - var_col_name, str: the column name for the variants
            ie 'V39D:D40G'
        - seq_col_name, str: the column name for the full sequence
        - fit_col_name, str: the column name for the fitness or yield
        - selectivity_col_name, str: the column name for the er of the enatiomer
        - protein_name, str: the protein name
        - seq_dir, str: the directory for the parent sequence fasta files
        - structure_dir, str: the directory for the structure files
        - mut_fasta_dir, str: the directory for the mutated fasta files
        """

        self._input_csv = input_csv
        self._scale_fit = scale_fit

        self._combo_col_name = combo_col_name
        self._var_col_name = var_col_name
        self._seq_col_name = seq_col_name
        self._fit_col_name = fit_col_name
        self._selectivity_col_name = selectivity_col_name
        self._protein_name = protein_name
        self._seq_dir = seq_dir
        self._structure_dir = structure_dir
        self._mut_fasta_dir = mut_fasta_dir

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
        if self._protein_name != "":
            return self._protein_name
        else:
            return self.lib_info["enzyme"]

    @property
    def parent_aa(self) -> str:
        """Return the parent amino acid"""
        return "".join(self.lib_info.get(self._combo_col_name, {}).values())

    @property
    def parent_fitness(self) -> float:
        """Return the parent fitness"""
        if self.parent_aa:
            parent_row = self.input_df[
                self.input_df[self._combo_col_name] == self.parent_aa
            ]
        parent_row = self.input_df[self.input_df[self._var_col_name] == "WT"]

        return parent_row[self._fit_col_name].values[0]

    @property
    def max_fitness(self) -> float:
        """Return the max fitness"""
        return self.input_df[self._fit_col_name].max()

    @property
    def norm_fit_factor(self) -> float:
        """Return the normalization factor"""
        if self._scale_fit == "parent":
            return self.parent_fitness
        elif self._scale_fit == "max":
            return self.max_fitness
        else:
            return 1

    @property
    def n_site(self) -> int:
        """Return the number of sites"""
        return len(self.lib_info["positions"])

    @property
    def mut_pos_list(self) -> list:
        """Return the list of positions"""
        return list(self.lib_info["positions"].values())

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
    def input_df_length(self):
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

    @property
    def structure_file(self) -> str:
        """Return the path to the structure file"""
        # find out if extention is cif or pdb
        if os.path.exists(os.path.join(self._structure_dir, f"{self.lib_name}.cif")):
            return os.path.join(self._structure_dir, f"{self.lib_name}.cif")
        else:
            return os.path.join(self._structure_dir, f"{self.lib_name}.pdb")

    @property
    def substrate_smiles(self) -> str:
        """Return the substrate smiles"""
        return canonicalize_smiles(self.lib_info["substrate-smiles"])

    @property
    def substrate_dets(self) -> str:
        """Return the substrate details"""
        return self.lib_info["substrate"]

    @property
    def cofactor_smiles(self) -> str:
        """Return the cofactor smiles"""
        return canonicalize_smiles(".".join(self.lib_info["cofactor-smiles"]))

    @property
    def cofactor_dets(self) -> str:
        """Return the cofactor details"""
        return "-".join(self.lib_info["cofactor"])


######### Handling SSM input meta data #########


class ProcessData(LibData):
    """
    A parent class to process the data
    """

    def __init__(
        self,
        input_csv: str,
        scale_fit: str = "",
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        selectivity_col_name: str = "er",
        protein_name: str = "",
        seq_dir: str = "data/seq",
        output_dir: str = "data/meta",
    ) -> None:

        """
        Args:
        - input_csv, str: path to the input csv file
        - scale_fit, str: ways to scale the fitness
            'parent' means the parent fitness = 1
            'max' means max fitness = 1
            '' or 'not_scaled' means no scaling, default
        - combo_col_name, str: the column name for the mutated amino acids
        - var_col_name, str: the column name for the variants
        - seq_col_name, str: the column name for the full sequence
        - fit_col_name, str: the column name for the fitness
        - selectivity_col_name, str: the column name for the er of the enatiomer
        - protein_name, str: the protein name
        - seq_dir, str: the directory for the parent sequence fasta files
        - output_dir, str: the output directory
        """

        super().__init__(
            input_csv=input_csv,
            scale_fit=scale_fit,
            combo_col_name=combo_col_name,
            var_col_name=var_col_name,
            seq_col_name=seq_col_name,
            fit_col_name=fit_col_name,
            selectivity_col_name=selectivity_col_name,
            protein_name=protein_name,
            seq_dir=seq_dir,
        )

        
        if self.scale_type == "not_scaled" or scale_fit == "":
            self._output_subdir = checkNgen_folder(
                os.path.join(output_dir)
            )
        else:
            self._output_subdir = checkNgen_folder(
                os.path.join(output_dir, self.scale_type)
            )

        print(f"Processing {self._input_csv} with {self._scale_fit}...")

        self._processed_df = self._process()

        # save the output csv
        self._processed_df.to_csv(self.output_csv, index=False)
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

        mut_name = ""
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
                    mut_name += ":" + wt_aa_pos + mut
                else:
                    mut_name += wt_aa_pos + mut

        return mut_name

    def _append_mut(self, df: pd.DataFrame) -> pd.DataFrame:

        """
        Apply the convert_muts function to the dataframe

        Args:
        - df, pd.DataFrame: the input dataframe

        Returns:
        - pd.DataFrame: the appended dataframe
        """

        df_appended = df.copy()
        df_appended.loc[:, self._var_col_name] = df_appended.apply(
            lambda x: self._convert_muts(x[self._combo_col_name]),
            axis=1,
        )

        df_appended[self._var_col_name] = df_appended[self._var_col_name].replace(
            "", "WT"
        )

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

    def _muts2combo(self, muts: str) -> str:
        """
        Convert the mutation with the form of parentaa1loc1mutaa1:parentaa2loc2mutaa2
        to the form of `VDGV` length equal to `n_site`
        with the missing position in self.mut_pos_list and
        corresponding wild type AA in self.parent_aa

        Args:
        - muts, str: the mutations that are 1 indexed
            ie. C1A:A3D
        
        Returns:
        - str: the mutated sequence length of mutated sites
        """

        if muts == "WT":
            return self.parent_aa

        combo_list = list(self.parent_aa)
        
        for mut in muts.split(":"):
            # first character is the parent amino acid and last is the mutated amino acid
            wt_aa, mut_aa = mut[0], mut[-1]

            # Extract position (1-indexed)
            mut_pos = int(mut[1:-1])

            # Get the index of the mutation in the parent_aa
            mut_pos_idx = self.mut_pos_list.index(mut_pos)

            # Assert that the parent sequence has the expected amino acid at the position
            assert (
                combo_list[mut_pos_idx] == wt_aa
            ), f"Mismatch at position {mut_pos}: expected {wt_aa}, found {combo_list[mut_pos_idx]}"

            # Replace the parent amino acid with the mutated amino acid at the specified position
            combo_list[mut_pos_idx] = mut_aa

        return "".join(combo_list)

    def _append_combo(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Apply the convert_muts function to the dataframe

        Args:
        - df, pd.DataFrame: the input dataframe

        Returns:
        - pd.DataFrame: the appended dataframe
        """

        df_appended = df.copy()
        df_appended.loc[:, self._combo_col_name] = df_appended.apply(
            lambda x: self._muts2combo(x[self._var_col_name]),
            axis=1,
        )

        return df_appended.copy()

    def _append_active_cutoff(self, df) -> pd.DataFrame:

        """
        Calculate the cutoff for active mutants based on
        1.96 standard deviations above the mean fitness of all stop-codon-containing sequences

        Returns:
        - pd.DataFrame: input dataframe with active column
        """

        if "active_cutoff" in self.lib_info:
            fit_cutoff = self.lib_info["active_cutoff"]

        elif df[self._var_col_name].str.contains("\*").any():

            stop_df = df[df[self._var_col_name].str.contains("\*")]
            avg_stop = stop_df[self._fit_col_name].mean()
            std_stop = stop_df[self._fit_col_name].std()
            fit_cutoff = 1.96 * std_stop + avg_stop

        else:
            return df

        # Apply the lambda functions to the DataFrame
        df["active"] = df[self._fit_col_name] > fit_cutoff

        return df

    def _process(self) -> pd.DataFrame:

        """
        Process the input csv
        """

        df = self.input_df.copy()

        # append muts column for none SSM data
        if self._var_col_name not in df.columns and self._combo_col_name in df.columns:
            df = self._append_mut(df).copy()

        # append combo column if only var data
        if self._var_col_name in df.columns and self._combo_col_name not in df.columns:
            df = self._append_combo(df).copy()

        # add mut number
        df["n_mut"] = df[self._var_col_name].str.split(":").str.len()

        # change WT n_mut to 0
        df.loc[df[self._var_col_name] == "WT", "n_mut"] = 0

        # add full seq from fasta file by modifying self.parent_seq with the mutations
        if self._seq_col_name not in df.columns and self._var_col_name in df.columns:
            df.loc[:, self._seq_col_name] = df[self._var_col_name].apply(
                lambda x: self._mut2seq(x)
            )

        # do more processing if not substrate scope data
        if "scope" not in self.lib_name:
            # scale the fitness
            df[self._fit_col_name] = df[self._fit_col_name] / self.norm_fit_factor


            # convert er to ee
            if ("er" in self._selectivity_col_name.lower()) and (
                self._selectivity_col_name in df.columns
            ):
                df["selectivity"] = df[self._selectivity_col_name].apply(er2ee)
                # drop the er column
                df.drop(self._selectivity_col_name, axis=1, inplace=True)


            # if there are multiple entries for the same combo name, take the mean
            groupby_cols = [col for col in df.columns if not pd.api.types.is_numeric_dtype(df[col])]
            
            df = df.groupby(groupby_cols).mean().reset_index()

            # split the amino acids for SSM data
            if "AA1" not in df.columns and self._combo_col_name in df.columns:
                df = self._split_aa(df).copy()

            # add active column
            df = self._append_active_cutoff(df).copy()

            # drop stop codon containing rows
            df = df[~df[self._var_col_name].str.contains("\*")].copy()

            # add col for enzyme name, substrate, cofactor, and their smile strings if relevant
            for col in APPEND_INFO_COLS:
                if col in self.lib_info:
                    append_info = self.lib_info[col]
                    # if the content is a list, convert it to a string
                    if isinstance(append_info, list):
                        df[col] = ".".join(append_info)
                    else:
                        df[col] = self.lib_info[col]

        # make n_mut to be int
        df["n_mut"] = df["n_mut"].astype(int)

        # make sure fitness is numeric
        if df[self._fit_col_name].dtype == "str":
            # convert fitness to numerical
            df[self._fit_col_name] = pd.to_numeric(df[self._fit_col_name].str.replace(",", ""), errors='coerce')

        # save the output csv
        return df

    @property
    def output_csv(self) -> str:

        """Return the path to the output csv"""

        return os.path.join(self._output_subdir, os.path.basename(self._input_csv))

    @property
    def processed_df(self) -> pd.DataFrame:

        """Return the output dataframe"""

        return self._processed_df


def preprocess_all(
    input_pattern: str = "data/lib/*.csv",
    scale_fit: str = "",
    combo_col_name: str = "AAs",
    var_col_name: str = "var",
    seq_col_name: str = "seq",
    fit_col_name: str = "fitness",
    selectivity_col_name: str = "er",
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
        '' or 'not_scaled' means no scaling, default
    - combo_col_name, str: the column name for the mutated amino acids
    - var_col_name, str: the column name for the variants
    - seq_col_name, str: the column name for the full sequence
    - fit_col_name, str: the column name for the fitness
    - seq_dir, str: the directory for the parent sequence fasta files
    - output_dir, str: the output directory

    """

    for input_csv in tqdm(sorted(glob(input_pattern))):

        ProcessData(
            input_csv=input_csv,
            scale_fit=scale_fit,
            combo_col_name=combo_col_name,
            var_col_name=var_col_name,
            seq_col_name=seq_col_name,
            fit_col_name=fit_col_name,
            selectivity_col_name=selectivity_col_name,
            seq_dir=seq_dir,
            output_dir=output_dir,
        )._process()


######### Handling ZS data #########


class ZSData(LibData):

    """
    A class for generating ZS scores
    """

    def __init__(
        self,
        input_csv: str,
        scale_fit: str = "",
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        mut_col_name: str = "mut",
        pos_col_name: str = "pos",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        sub_col_name: str = "substrate",
        subsmiles_col_name: str = "substrate_smiles",
        cof_col_name: str = "cofactor",
        cofsmiles_col_name: str = "cofactor_smiles",
        joint_col_name: str = "joint",
        jointsmiles_col_name: str = "joint_smiles",
        protein_name: str = "",
        withsub: bool = True,
        seq_dir: str = "data/seq",
        structure_dir: str = "data/structure",
        zs_dir: str = "zs",
    ):

        """
        - mut_col_name, str: the column name for the mutations
            ie ['A', 'D']
        - pos_col_name, str: the column name for the positions
            ie [39, 40]

        input_csv: str,
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        protein_name: str = "",
        seq_dir: str = "data/seq",
        structure_dir: str = "data/structure",
        mut_fasta_dir: str = "data/mut_fasta",
        """

        super().__init__(
            input_csv=input_csv,
            scale_fit=scale_fit,
            combo_col_name=combo_col_name,
            var_col_name=var_col_name,
            seq_col_name=seq_col_name,
            fit_col_name=fit_col_name,
            protein_name=protein_name,
            seq_dir=seq_dir,
            structure_dir=structure_dir,
        )

        self._mut_col_name = mut_col_name
        self._pos_col_name = pos_col_name
        self._sub_col_name = sub_col_name
        self._subsmiles_col_name = subsmiles_col_name
        self._cof_col_name = cof_col_name
        self._cofsmiles_col_name = cofsmiles_col_name
        self._joint_col_name = joint_col_name
        self._withsub = withsub

        self._zs_dir = checkNgen_folder(zs_dir)

    def _append_mut_dets(self, var: str) -> tuple:

        """
        Append mut details from the combo column

        Args:
        - var, str: the variants sequence

        Returns:
        - list: the list of mutated AA
        - list: the list of mutated positions
        """

        mut_list = []
        pos_list = []

        if var == "WT":
            return mut_list, pos_list

        for v in var.split(":"):

            mut_list.append(v[-1])
            # note the info dict positiosn is 1 indexed
            pos_list.append(int(v[1:-1]))

        return mut_list, pos_list

    @property
    def df(self) -> pd.DataFrame:

        """
        Get the dataframe with mutation details
        """

        df = self.input_df.copy()

        df[[self._mut_col_name, self._pos_col_name]] = df.apply(
            lambda x: pd.Series(self._append_mut_dets(x[self._var_col_name])),
            axis=1,
        )

        # filter out stop codon if any

        return df[~df[self._var_col_name].str.contains("\*")].copy()

    @property
    def max_n_mut(self) -> int:

        """
        Get the maximum number of mutations
        """

        return self.df["n_mut"].max()

    @property
    def zs_struct_name(self) -> str:

        """
        Get the name of the structure file
        """

        if self._withsub:
            return self.lib_name
        else:
            return self.lib_info["enzyme"]


def gen_apo_structures(
    struct_dir: str = "data/structure",
    chain_id: str = "A",
    residues_to_remove: list = ["HOH", "SO4"],
) -> None:

    """
    Generate cleaned up (no water, SO4) and apo structures for all the protein structures taken from PDB
    and save to the apo subdirectory
    """

    clean_dir = checkNgen_folder(os.path.join(struct_dir, "clean"))
    apo_dir = checkNgen_folder(os.path.join(struct_dir, "apo"))

    # get all the pdb files in the structure directory
    pdb_files = sorted(glob(os.path.join(struct_dir, "*.pdb")))

    for input_struct in pdb_files:

        protein_name = get_file_name(input_struct)

        clean_struct = os.path.join(clean_dir, f"{protein_name}.pdb")
        apo_struct = os.path.join(apo_dir, f"{protein_name}.pdb")
        temp_path = apo_struct.replace(".pdb", "_temp.pdb")

        # get the name of the structure file
        get_chain_structure(
            input_file_path=input_struct,
            output_file_path=temp_path,
            chain_id=chain_id,
        )

        remove_residues_from_pdb(
            input_pdb=input_struct,
            output_pdb=clean_struct,
            residues_to_remove=residues_to_remove,
        )

        remove_hetatm(
            input_pdb=temp_path,
            output_pdb=apo_struct,
        )

        os.remove(temp_path)


def process_substratescope_df(
    pattern: str = "data/substrate_scope/*",
    output_dir: str = "data/meta",
    var_col_name: str = "var",
    ):
    """
    Process the substrate scope data
    """

    output_dir = checkNgen_folder(output_dir)

    for input_csv in tqdm(sorted(glob(pattern))):

        df = pd.read_csv(input_csv)

        df["n_mut"] = df[var_col_name].str.split(":").str.len()

        # change WT n_mut to 0
        df.loc[df[var_col_name] == "WT", "n_mut"] = 0

        # save the dataframe to a csv file
        output_csv = os.path.join(output_dir, f"{get_file_name(input_csv)}.csv")
        df.to_csv(output_csv, index=False)