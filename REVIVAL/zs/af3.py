"""
A script to generate JSON files for AlphaFold3 inference from a CSV file.
"""

from __future__ import annotations
import concurrent.futures

import os
import json
import subprocess

from copy import deepcopy
from glob import glob

from tqdm import tqdm

from REVIVAL.preprocess import ZSData
from REVIVAL.util import checkNgen_folder, canonicalize_smiles, load_json


class AF3Prep(ZSData):
    def __init__(
        self,
        input_csv: str,
        scale_fit: str = "not_scaled",
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        mut_col_name: str = "mut",
        pos_col_name: str = "pos",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        seq_dir: str = "data/seq",
        zs_dir: str = "zs",
        af3_dir: str = "af3",
        gen_opt: str = "joint",
        cofactor_dets: str = "cofactor",
        ifrerun: bool = False,
        full_model_path: str = "/disk2/fli/af3_inference/model",
        full_db_path: str = "/disk2/fli/af3_inference/databases",
        max_workers: int = 16,  # for msa generation parallelization
        gpu_id: str = "0",
    ):

        super().__init__(
            input_csv=input_csv,
            scale_fit=scale_fit,
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

        self._ifrerun = ifrerun

        self._full_model_path = full_model_path
        self._full_db_path = full_db_path

        self._gpu_id = gpu_id

        self._af3_dir = checkNgen_folder(os.path.join(self._zs_dir, f"{af3_dir}"))

        self._af3_msa_insubdir = checkNgen_folder(
            os.path.join(self._af3_dir, "msa_in", self.protein_name)
        )
        self._af3_msa_subdir = checkNgen_folder(
            os.path.join(self._af3_dir, "msa", self.protein_name)
        )

        self._af3_json_subdir = checkNgen_folder(
            os.path.join(self._af3_dir, f"json_{self._gen_opt}", self.lib_name)
        )

        self._af3_struct_subdir = checkNgen_folder(
            os.path.join(self._af3_dir, f"struct_{self._gen_opt}", self.lib_name)
        )

        self._sub_smiles = canonicalize_smiles(self.lib_info["substrate-smiles"])
        self._sub_dets = self.lib_info["substrate"]

        self._cofactor_smiles = canonicalize_smiles(
            ".".join(self.lib_info[f"{cofactor_dets}-smiles"])
        )
        self._cofactor_dets = "-".join(self.lib_info[cofactor_dets])

        self._joint_smiles = self._sub_smiles + "." + self._cofactor_smiles
        self._joint_dets = self._sub_dets + "_" + self._cofactor_dets

        # Parallelizing MSA generation only
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=max_workers
        ) as executor:
            futures = []
            for i, (var, seq) in enumerate(
                self.df[[self._var_col_name, self._seq_col_name]].values
            ):
                futures.append(
                    executor.submit(
                        self._gen_var_msa, var=var.replace(":", "_"), seq=seq
                    )
                )

            for future in concurrent.futures.as_completed(futures):
                try:
                    future.result()
                except Exception as e:
                    print(f"Error in parallel execution: {e}")

        # now do the inferences and the _gen_var_msa step will simply check the path
        for (
            var,
            seq,
        ) in tqdm(self.df[[self._var_col_name, self._seq_col_name]].values):

            rename_var = var.replace(":", "_")

            var_msa_outpath = self._gen_var_msa(var=rename_var, seq=seq)
            self._run_var_inference(
                var=rename_var, seq=seq, var_msa_outpath=var_msa_outpath
            )

    def _gen_var_msa(self, var, seq):

        """
        A method to generate the MSA for each var, seq,
        which can then be used for different variants.
        """

        var_msa_inpath = os.path.join(self._af3_msa_insubdir, f"{var}.json")

        # all cases will be lower and all : will be stripped
        modified_var = var.lower().replace(":", "")

        # _data is added by af3
        var_msa_outpath = os.path.join(
            self._af3_msa_subdir, modified_var, f"{modified_var}_data.json"
        )

        # check if exist and ifrerun is False
        if os.path.exists(var_msa_outpath) and not self._ifrerun:
            print(f"MSA for {var} already exists at {var_msa_outpath}. Skipping...")
            return var_msa_outpath

        json_data = {
            "name": f"{var}",
            "sequences": [{"protein": {"id": "A", "sequence": f"{seq}"}}],
            "modelSeeds": [1],
            "dialect": "alphafold3",
            "version": 1,
        }

        # Save the JSON to generate msa
        with open(var_msa_inpath, "w") as json_file:
            json.dump(json_data, json_file, indent=4)

        # run the docker command to generate the msa
        docker_command = [
            "sudo",
            "docker",
            "run",
            # "-it",
            # "--env", f"CUDA_VISIBLE_DEVICES={self._gpu_id}",
            "--volume",
            f"{os.path.abspath(os.path.dirname(var_msa_inpath))}:/root/af_input",
            "--volume",
            f"{os.path.abspath(os.path.dirname(os.path.dirname(var_msa_outpath)))}:/root/af_output",
            "--volume",
            f"{self._full_model_path}:/root/models",
            "--volume",
            f"{self._full_db_path}:/root/public_databases",
            "--gpus",
            f"device={self._gpu_id}",
            # "all",
            "alphafold3",
            "python",
            "run_alphafold.py",
            f"--json_path=/root/af_input/{os.path.basename(var_msa_inpath)}",
            "--jackhmmer_n_cpu=8",
            "--nhmmer_n_cpu=8",
            "--model_dir=/root/models",
            "--output_dir=/root/af_output",
            "--run_data_pipeline=True",
            "--run_inference=False",
        ]

        # Log file for the output
        log_file_path = var_msa_inpath.replace(".json", ".log")

        try:
            # Open the log file
            with open(log_file_path, "w") as log_file:
                # Run the Docker command using subprocess
                result = subprocess.run(
                    docker_command,
                    check=True,
                    text=True,
                    stdout=log_file,  # Redirect stdout to log file
                    stderr=log_file,  # Redirect stderr to log file
                )
            print(
                f"Docker command executed successfully. Logs saved to {log_file_path}"
            )
        except subprocess.CalledProcessError as e:
            print(
                f"Error occurred while executing Docker command. See logs in {log_file_path}"
            )
            print(f"Error: {e}")

        return var_msa_outpath

    def _run_var_inference(self, var, seq, var_msa_outpath):

        """
        A method to generate the AlphaFold3 JSON files for each variant.
        """

        var_dir = os.path.join(self._af3_struct_subdir, var.lower())

        # the output cif file
        model_cif = os.path.join(var_dir, f"{var.lower()}_model.cif")

        # del the folder if the model is not there but the folder exists
        if os.path.exists(var_dir) and not os.path.exists((model_cif)):
            os.system(f"sudo rm -r {var_dir}")

        # check if the output model exists and ifrerun is False
        if (
            os.path.exists(model_cif)
            and not self._ifrerun
        ):
            print(
                f"Structures for {var} already exist at {self._af3_struct_subdir}. Skipping..."
            )
            return

        msa_json = load_json(var_msa_outpath).get("sequences", [])[0].get("protein", {})

        # Create the JSON structure
        json_data = {
            "name": f"{var}",
            "sequences": [
                {
                    "protein": {
                        "id": "A",
                        "sequence": f"{seq}",
                        "modifications": msa_json.get("modifications", None),
                        "unpairedMsa": msa_json.get("unpairedMsa", None),
                        "pairedMsa": msa_json.get("pairedMsa", None),
                        "templates": msa_json.get("templates", None),
                    }
                }
            ],
            "modelSeeds": [1],
            "dialect": "alphafold3",
            "version": 2,
        }

        if self._gen_opt == "substrate-no-cofactor":

            # add substrate
            json_data["sequences"].append(
                {"ligand": {"id": "B", "smiles": f"{self._sub_smiles}"}}
            )

        elif self._gen_opt == "joint-cofactor-no-substrate":

            json_data["sequences"].append(
                {"ligand": {"id": "C", "smiles": f"{self._cofactor_smiles}"}}
            )

        elif self._gen_opt == "joint-cofactor-seperate-substrate":

            # add substrate
            json_data["sequences"].append(
                {"ligand": {"id": "B", "smiles": f"{self._sub_smiles}"}}
            )

            # add cofactor
            json_data["sequences"].append(
                {"ligand": {"id": "C", "smiles": f"{self._cofactor_smiles}"}}
            )

        elif self._gen_opt == "seperate":

            # add substrate
            json_data["sequences"].append(
                {"ligand": {"id": "B", "smiles": f"{self._sub_smiles}"}}
            )

            for j, (cofactor_dets, cofactor_smiles) in enumerate(
                zip(self.lib_info["cofactor"], self.lib_info["cofactor-smiles"])
            ):
                cofactor_smiles = canonicalize_smiles(cofactor_smiles)
                json_data["sequences"].append(
                    {
                        "ligand": {
                            "id": chr(
                                67 + j
                            ),  # 'C', 'D', etc. for additional cofactors
                            "smiles": cofactor_smiles,
                        }
                    }
                )

        else:

            json_data["sequences"].append(
                {"ligand": {"id": "B", "smiles": f"{self._joint_smiles}"}}
            )

        # Save each JSON to a separate file
        json_file_path = os.path.join(self._af3_json_subdir, f"{var}.json")
        with open(json_file_path, "w") as json_file:
            json.dump(json_data, json_file, indent=4)
            print(f"JSON file saved to {json_file_path}")

        # run the docker command
        docker_command = [
            "sudo",
            "docker",
            "run",
            # "-it",
            # "--env", f"CUDA_VISIBLE_DEVICES={self._gpu_id}",
            "--volume",
            f"{os.path.abspath(os.path.dirname(json_file_path))}:/root/af_input",
            "--volume",
            f"{os.path.abspath(self._af3_struct_subdir)}:/root/af_output",
            "--volume",
            f"{self._full_model_path}:/root/models",
            "--volume",
            f"{self._full_db_path}:/root/public_databases",
            "--gpus",
            f"device={self._gpu_id}",
            # "all",
            "alphafold3",
            "python",
            "run_alphafold.py",
            f"--json_path=/root/af_input/{os.path.basename(json_file_path)}",
            "--model_dir=/root/models",
            "--output_dir=/root/af_output",
            "--run_data_pipeline=False",
            "--run_inference=True",
        ]

        # Log file for the output
        log_file_path = json_file_path.replace(".json", ".log")

        try:
            # Open the log file
            with open(log_file_path, "w") as log_file:
                # Run the Docker command using subprocess
                result = subprocess.run(
                    docker_command,
                    check=True,
                    text=True,
                    stdout=log_file,  # Redirect stdout to log file
                    stderr=log_file,  # Redirect stderr to log file
                )
            print(
                f"Docker command executed successfully. Logs saved to {log_file_path}"
            )
        except subprocess.CalledProcessError as e:
            print(
                f"Error occurred while executing Docker command. See logs in {log_file_path}"
            )
            print(f"Error: {e}")


def run_af3_prep(
    pattern: str | list = "data/meta/not_scaled/*.csv",
    gen_opt: str = "joint",
    cofactor_dets: str = "cofactor",
    gpu_id: str = "0",
    kwargs: dict = {},
):
    """
    Run af3 for all libraries

    Args:
    - pattern: str | list: the pattern for the input csv files
    - gen_opt: str: The generation option for the AF3 structure
    - kwargs: dict: The arguments for the AF3Data class
    """

    if isinstance(pattern, str):
        lib_list = sorted(glob(pattern))
    else:
        lib_list = deepcopy(pattern)

    for lib in lib_list:
        print(f"Running AF3 for {lib}...")
        AF3Prep(input_csv=lib, gen_opt=gen_opt, cofactor_dets=cofactor_dets, **kwargs)