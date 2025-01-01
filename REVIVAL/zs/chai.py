"""
A script for chai structures and also extract the score.
"""

from __future__ import annotations

import os
from glob import glob
from copy import deepcopy
from tqdm import tqdm
from pathlib import Path

import numpy as np
import pandas as pd

import torch

from rdkit import Chem

from chai_lab.chai1 import run_inference

from REVIVAL.preprocess import ZSData
from REVIVAL.util import checkNgen_folder, get_file_name


def canonicalize_smiles(smiles_string: str) -> str:

    """
    A function to canonicalize a SMILES string.

    Args:
    - smiles_string (str): The input SMILES string.

    Returns:
    - str: The canonicalized SMILES string.
    """

    molecule = Chem.MolFromSmiles(smiles_string)
    if molecule:
        canonical_smiles = Chem.MolToSmiles(molecule, canonical=True)
        return canonical_smiles


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
        gen_opt: str = "joint",
        ifrerun: bool = False,
        torch_device: str = "cuda",
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

        self._gen_opt = gen_opt

        self._chai_dir = checkNgen_folder(os.path.join(self._zs_dir, f"{chai_dir}_{self._gen_opt}"))
        self._chai_struct_dir = checkNgen_folder(
            os.path.join(self._chai_dir, chai_struct_dir)
        )
        self._chai_struct_subdir = checkNgen_folder(
            os.path.join(self._chai_struct_dir, self.lib_name)
        )

        self._ifrerun = ifrerun
        self._torch_device = torch_device

        self._gen_chai_structure()

    def _gen_chai_structure(self):
        """
        A method to generate the chai structure for each variant.
        """

        for (
            var,
            seq,
        ) in tqdm(self.df[[self._var_col_name, self._seq_col_name]].values):

            output_subdir = os.path.join(self._chai_struct_subdir, var)

            # Need to clean up the sequence
            seq = seq.strip().replace("*", "").replace(" ", "").upper()

            input_fasta = f">protein|{self.lib_name}_{var}\n{seq}\n"

            sub_smiles = canonicalize_smiles(self.lib_info["substrate-smiles"])
            sub_dets = self.lib_info["substrate"]

            cofactor_smiles = canonicalize_smiles(".".join(self.lib_info["cofactor-smiles"]))
            cofactor_dets = "-".join(self.lib_info["cofactor"])

            joint_smiles = sub_smiles + "." + cofactor_smiles
            joint_dets = sub_dets + "_" + cofactor_dets

            if self._gen_opt == "joint-cofactor-no-substrate":

                # now add cofactor
                input_fasta += f">ligand|{cofactor_dets}\n{cofactor_smiles}\n"

            elif self._gen_opt == "joint-cofactor-seperate-substrate":

                # add substrate first
                input_fasta += f">ligand|{sub_dets}\n{sub_smiles}\n"

                # now get cofactor smiles
                input_fasta += f">ligand|{cofactor_dets}\n{cofactor_smiles}\n"

            elif self._gen_opt == "substrate-no-cofactor":

                # now add substrate
                input_fasta += f">ligand|{sub_dets}\n{sub_smiles}\n"

            else:
            
                # now add substrate
                input_fasta += f">ligand|{joint_smiles}\n{joint_dets}\n"

            # only rerun if the flag is set and the output folder doies not exists
            if self._ifrerun or not os.path.exists(output_subdir):

                output_subdir = Path(checkNgen_folder(output_subdir))

                fasta_path = Path(f"{output_subdir}/{var}.fasta")
                fasta_path.write_text(input_fasta)

                print(f"Running chai for {var}...")

                run_inference(
                    fasta_file=fasta_path,
                    output_dir=output_subdir,
                    # 'default' setup
                    num_trunk_recycles=3,
                    num_diffn_timesteps=200,
                    seed=42,
                    device=torch.device(self._torch_device),
                    use_esm_embeddings=True,
                )

                renamed_output_files = []

                # get name of the output cif or pdb files
                output_strcut_files = sorted(
                    glob(f"{output_subdir}/*.cif") + glob(f"{output_subdir}/*.pdb")
                )

                # rename the output files cif or pdb files
                for output_strcut_file in output_strcut_files:
                    renamed_output_file = output_strcut_file.replace(
                        "pred.model_idx", var
                    )
                    os.rename(
                        output_strcut_file,
                        renamed_output_file,
                    )
                    renamed_output_files.append(renamed_output_file)

                renamed_scores_files = []

                # for npz files do the same
                output_scores_files = sorted(glob(f"{output_subdir}/*.npz"))

                for output_scores_file in output_scores_files:
                    renamed_output_file = output_scores_file.replace(
                        "scores.model_idx", var
                    )
                    os.rename(
                        output_scores_file,
                        renamed_output_file,
                    )
                    renamed_scores_files.append(renamed_output_file)

            else:
                print(f"{var} exists and ifrerun {self._ifrerun}...")
                renamed_output_files = glob(f"{output_subdir}/*.cif") + glob(
                    f"{output_subdir}/*.pdb"
                )
                renamed_scores_files = glob(f"{output_subdir}/*.npz")

        return renamed_output_files, renamed_scores_files


def run_gen_chai_structure(
    pattern: str | list = "data/meta/scale2parent/*.csv", 
    gen_opt: str = "joint",
    kwargs: dict = {}
):
    """
    Run the chai gen mut file function for all libraries

    Args:
    - pattern: str | list: the pattern for the input csv files
    - gen_opt: str: The generation option for the chai structure
    - kwargs: dict: The arguments for the ChaiData class
    """

    if isinstance(pattern, str):
        lib_list = glob(pattern)
    else:
        lib_list = deepcopy(pattern)

    for lib in lib_list:
        print(f"Running chai gen mut file for {lib}...")
        ChaiData(input_csv=lib, gen_opt=gen_opt, **kwargs)


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

    for subfolder in sorted(glob(f"{mut_structure_dir}/*")):
        var = os.path.basename(subfolder)

        for rep_npz in sorted(glob(f"{subfolder}/*.npz")):

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