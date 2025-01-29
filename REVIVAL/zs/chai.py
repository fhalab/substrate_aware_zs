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

from Bio.PDB import MMCIFParser

import torch

from chai_lab.chai1 import run_inference

from REVIVAL.preprocess import ZSData
from REVIVAL.global_param import LIB_INFO_DICT
from REVIVAL.chem_helper import canonicalize_smiles
from REVIVAL.util import checkNgen_folder


def get_residue_uncertainty(cif_file, chain_id, resid_list):
    """
    Extract and compute the average B-factors (uncertainties) for atoms in a given list of residues.

    Args:
        cif_file (str): Path to the .cif file.
        chain_id (str): Chain ID of the residues.
        resid_list (list): List of residue IDs.

    Returns:
        dict: Residue IDs with their average B-factors (uncertainties).
        float: Overall average B-factor across all residues.
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("protein", cif_file)

    residue_bfactors = {}  # To store average B-factors for each residue
    all_bfactors = []  # To calculate the overall average

    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    if residue.id[1] in resid_list:
                        bfactors = [atom.bfactor for atom in residue]
                        avg_bfactor = sum(bfactors) / len(bfactors)
                        residue_bfactors[residue.id[1]] = avg_bfactor
                        all_bfactors.extend(bfactors)

    if not residue_bfactors:
        raise ValueError(
            f"No residues in {resid_list} found in chain {chain_id} of {cif_file}."
        )

    return sum(all_bfactors) / len(all_bfactors)


class ChaiStruct(ZSData):
    def __init__(
        self,
        input_csv: str,
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        mut_col_name: str = "mut",
        pos_col_name: str = "pos",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        seq_dir: str = "data/seq",
        zs_dir: str = "zs",
        chai_dir: str = "chai",
        chai_struct_dir: str = "struct",
        gen_opt: str = "joint",
        cofactor_dets: str = "cofactor",
        ifrerun: bool = False,
        samesub: bool = True,
        torch_device: str = "cuda",
    ):

        """
        A class to generate the chai structure for each variant.

        Args:
        - input_csv, str: The path to the input csv file
        - combo_col_name, str: The column name for the combo
        - var_col_name, str: The column name for the variant
        - mut_col_name, str: The column name for the mutation
        - pos_col_name, str: The column name for the position
        - seq_col_name, str: The column name for the sequence
        - fit_col_name, str: The column name for the fitness
        - seq_dir, str: The path to the sequence directory
        - zs_dir, str: The path to the zs directory
        - chai_dir, str: The path to the chai directory
        - chai_struct_dir, str: The path to the chai structure directory
        - gen_opt, str: The generation option for the chai structure
        - cofactor_dets, str: The cofactor details, ie "cofactor" or "inactivated-cofactor"
        - ifrerun, bool: A flag to rerun the chai structure
        - torch_device, str: The torch device to use
        """

        super().__init__(
            input_csv=input_csv,
            combo_col_name=combo_col_name,
            var_col_name=var_col_name,
            mut_col_name=mut_col_name,
            pos_col_name=pos_col_name,
            seq_col_name=seq_col_name,
            fit_col_name=fit_col_name,
            seq_dir=seq_dir,
            zs_dir=zs_dir,
        )

        self._gen_opt = gen_opt

        self._samesub = samesub
        if "substratescope" in input_csv:
            self._samesub = False

        self._chai_dir = checkNgen_folder(os.path.join(self._zs_dir, chai_dir))
        self._chai_struct_dir = checkNgen_folder(
            os.path.join(self._chai_dir, f"{chai_struct_dir}_{self._gen_opt}")
        )

        append_dets = ""
        if cofactor_dets != "cofactor":
            append_dets = f"_{cofactor_dets}"

        self._chai_struct_subdir = checkNgen_folder(
            os.path.join(self._chai_struct_dir, f"{self.lib_name}{append_dets}")
        )

        self._ifrerun = ifrerun
        self._torch_device = torch_device

        # TODO can make the following more concise
        if samesub:

            self._cofactor_smiles = canonicalize_smiles(
                ".".join(self.lib_info[f"{cofactor_dets}-smiles"])
            )
            self._cofactor_dets = "-".join(self.lib_info[cofactor_dets])

            self._joint_smiles = self.substrate_smiles + "." + self._cofactor_smiles
            self._joint_dets = self.substrate_dets + "_" + self._cofactor_dets
            self._gen_chai_structure()
        else:
            self._gen_chai_structure_diff_sub()

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

            if self._gen_opt in ["no-substrate-no-cofactor", "apo", "empty"]:

                # do nothing
                pass

            elif self._gen_opt == "substrate-no-cofactor":

                # add substrate
                input_fasta += f">ligand|{self.substrate_dets}\n{self.substrate_smiles}\n"

            elif self._gen_opt == "joint-cofactor-no-substrate":

                # add cofactor
                input_fasta += (
                    f">ligand|{self._cofactor_dets}\n{self._cofactor_smiles}\n"
                )

            elif self._gen_opt == "joint-cofactor-seperate-substrate":

                # add substrate first
                input_fasta += f">ligand|{self.substrate_dets}\n{self.substrate_smiles}\n"

                # add cofactor smiles
                input_fasta += (
                    f">ligand|{self._cofactor_dets}\n{self._cofactor_smiles}\n"
                )

            elif self._gen_opt == "seperate":

                # add joint
                input_fasta += f">ligand|{self.substrate_dets}\n{self.substrate_smiles}\n"

                # loop through the cofactors and add them individually
                for cofactor_dets, cofactor_smiles in zip(
                    self.lib_info["cofactor"], self.lib_info["cofactor-smiles"]
                ):
                    cofactor_smiles = canonicalize_smiles(cofactor_smiles)
                    input_fasta += f">ligand|{cofactor_dets}\n{cofactor_smiles}\n"

            else:

                # add substrate
                input_fasta += f">ligand|{self._joint_dets}\n{self._joint_smiles}\n"

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

    def _gen_chai_structure_diff_sub(self):
        """
        A method to generate the chai structure for each variant.

        Assume additional columns: substrate, substrate-smiles, cofactor, cofactor-smiles, rxn_id
        """

        for (var, seq, sub, sub_smile, cof, cof_smile, rxn_id) in tqdm(
            self.df[
                [
                    self._var_col_name,
                    self._seq_col_name,
                    "substrate",
                    "substrate-smiles",
                    "cofactor",
                    "cofactor-smiles",
                    "rxn_id",
                ]
            ].values
        ):

            output_subdir = os.path.join(
                self._chai_struct_subdir, f"{var}_{str(rxn_id)}"
            )

            # Need to clean up the sequence
            seq = seq.strip().replace("*", "").replace(" ", "").upper()

            input_fasta = f">protein|{self.lib_name}_{var}_{rxn_id}\n{seq}\n"

            if self._gen_opt in ["no-substrate-no-cofactor", "apo", "empty"]:

                # do nothing
                pass

            elif self._gen_opt == "substrate-no-cofactor":

                # add substrate
                input_fasta += f">ligand|{sub}\n{sub_smile}\n"

            elif self._gen_opt == "joint-cofactor-no-substrate":

                # add cofactor
                input_fasta += f">ligand|{cof}\n{cof_smile}\n"

            elif self._gen_opt == "seperate":

                # add substrate first
                input_fasta += f">ligand|{sub}\n{sub_smile}\n"

                # add cofactor smiles
                input_fasta += f">ligand|{cof}\n{cof_smile}\n"

            else:
                joint_dets = sub + "_" + cof
                joint_smiles = sub_smile + "." + cof_smile

                # add substrate
                input_fasta += f">ligand|{joint_dets}\n{joint_smiles}\n"

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
    pattern: str | list = "data/meta/not_scaled/*.csv",
    gen_opt: str = "joint",
    cofactor_dets: str = "cofactor",
    samesub: bool = True,
    kwargs: dict = {},
):
    """
    Run the chai gen mut file function for all libraries

    Args:
    - pattern: str | list: the pattern for the input csv files
    - gen_opt: str: The generation option for the chai structure
    - cofactor_dets: str: The cofactor details, ie "cofactor" or "inactivated-cofactor"
    - kwargs: dict: The arguments for the ChaiStruct class
    """

    if isinstance(pattern, str):
        lib_list = glob(pattern)
    else:
        lib_list = deepcopy(pattern)

    for lib in lib_list:
        print(f"Running chai gen mut file for {lib}...")
        ChaiStruct(
            input_csv=lib,
            gen_opt=gen_opt,
            cofactor_dets=cofactor_dets,
            samesub=samesub,
            **kwargs,
        )


def parse_chai_scores(mut_structure_dir: str, score_dir_name: str = "score"):

    """
    A function for going through the subfolder and getting the chai scores
    to generate a dataframe with the following columns:
        - var: The mutation, ie I165A:I183A:Y301V
        - rep: The replicate number
        - aggregate_score
        - ptm
        - iptm
        - has_inter_chain_clashes

    Args:
    - input_dir, str: The path to the folder containing the chai score
        ie zs/chai/struct_joint/PfTrpB-4bromo
    - score_dir_name, str: The name of the score directory but keep other details
        ie zs/chai/score_joint
    """

    # Determine the docking details based on folder naming
    dock_dets = (
        "joint"
        if "joint" in mut_structure_dir
        else "separate"
        if "seperate" in mut_structure_dir
        else ""
    )

    output_dir = checkNgen_folder(
        os.path.dirname(mut_structure_dir).replace("struct", score_dir_name)
    )
    lib_name = os.path.basename(mut_structure_dir)

    # Prepare the score keys for data extraction
    overall_keys = ["aggregate_score", "ptm", "iptm"]
    chain_labels = ["A", "B", "C"] if dock_dets == "separate" else ["A", "B"]
    score_keys = deepcopy(overall_keys)
    for chain in chain_labels:
        score_keys.append(f"chain_ptm_{chain}")
        for other_chain in chain_labels:
            if chain != other_chain:
                score_keys.append(f"chain_iptm_{chain}{other_chain}")

    # Initialize results storage
    results = []

    # Iterate through each variant directory
    for subfolder in sorted(glob(os.path.join(mut_structure_dir, "*"))):
        var_name = os.path.basename(subfolder)
        var_data = {"var": var_name}

        # Prepare for averaging
        score_sums = {key: [] for key in score_keys}

        # Loop through each replicate score file
        for rep_index, rep_npz in enumerate(
            sorted(glob(os.path.join(subfolder, "*.npz")))
        ):
            try:
                npz = np.load(rep_npz)
                # Extract the score data for this replicate
                for key in overall_keys:
                    var_data[f"{key}_{rep_index}"] = npz[key][0]

                # add mean_site_score
                var_data[f"mean_site_score_{rep_index}"] = get_residue_uncertainty(
                    cif_file=rep_npz.replace(".npz", ".cif"),
                    chain_id="A",
                    resid_list=list(
                        LIB_INFO_DICT[os.path.basename(mut_structure_dir)][
                            "positions"
                        ].values()
                    ),
                )

                # Process chain-level ptm and iptm scores
                for i, chain in enumerate(chain_labels):
                    var_data[f"chain_ptm_{chain}_{rep_index}"] = npz["per_chain_ptm"][
                        0
                    ][i]
                    for j, other_chain in enumerate(chain_labels):
                        if i != j:
                            var_data[
                                f"chain_iptm_{chain}{other_chain}_{rep_index}"
                            ] = npz["per_chain_pair_iptm"][0][i, j]

                # Only consider for averaging if no clashes
                if not npz["has_inter_chain_clashes"][0]:
                    for key in score_keys:
                        value = var_data.get(f"{key}_{rep_index}")
                        if value is not None:
                            score_sums[key].append(value)

            except Exception as e:
                print(f"Error processing {rep_npz}: {e}")

        # Compute averages for the variant and store them as columns
        for key, values in score_sums.items():
            var_data[f"{key}_avg"] = np.mean(values) if values else None
            var_data[f"{key}_std"] = np.std(values) if values else None

        # add mean site score
        mean_site_scores = [
            var_data[f"mean_site_score_{rep_index}"]
            for rep_index in range(5)
            if var_data.get(f"mean_site_score_{rep_index}") is not None
        ]
        var_data["mean_site_score_avg"] = np.mean(mean_site_scores)
        var_data["mean_site_score_std"] = np.std(mean_site_scores)

        # Collect results
        results.append(var_data)

    # Convert results to a DataFrame and save it
    df = pd.DataFrame(results)

    df.to_csv(f"{output_dir}/{lib_name}.csv", index=False)
    print(f"Saved chai scores for {lib_name} to {output_dir}/{lib_name}.csv")


def parse_all_chai_scores(chai_struct_dir: str = "zs/chai/struct_joint"):
    """
    A function to parse all the chai scores for all libraries
    """

    for lib in tqdm(sorted(glob(f"{chai_struct_dir}/*"))):
        print(f"Parsing chai scores for {lib}...")
        parse_chai_scores(lib)