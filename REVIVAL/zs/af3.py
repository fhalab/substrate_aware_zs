"""
A script to generate JSON files for AlphaFold3 inference from a CSV file.
"""

from __future__ import annotations

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
        # af3_msa_indir: str = "msa_in",
        # af3_msa_dir: str = "msa",
        # af3_json_dir: str = "json",
        # af3_command_dir: str = "command",
        # af3_struct_dir: str = "struct",
        gen_opt: str = "joint",
        cofactor_dets: str = "cofactor",
        ifrerun: bool = False,
        full_model_path: str = "/disk2/fli/af3_inference/model",
        full_db_path: str = "/disk2/fli/af3_inference/databases",
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

        # self._af3_json_dir = af3_json_dir
        # self._af3_command_dir = af3_command_dir
        # self._af3_struct_dir = af3_struct_dir

        self._af3_dir = checkNgen_folder(
            os.path.join(self._zs_dir, f"{af3_dir}_{self._gen_opt}")
        )

        self._af3_msa_indir = checkNgen_folder(os.path.join(self._af3_dir, "msa_in"))
        self._af3_msa_inpath = os.path.join(
            self._af3_msa_indir, f"{self.protein_name}.json"
        )
        print(self._af3_msa_inpath, self.protein_name)
        self._af3_msa_exe = os.path.join(self._af3_msa_indir, f"{self.protein_name}.sh")

        self._af3_msa_dir = checkNgen_folder(os.path.join(self._af3_dir, "msa"))
        # _data is added by af3
        self._af3_msa_outpath = os.path.join(
            self._af3_msa_dir, f"{self.protein_name}_data.json"
        )

        self._af3_json_subdir = checkNgen_folder(
            os.path.join(self._af3_dir, "json", self.lib_name)
        )

        self._af3_command_subdir = checkNgen_folder(
            os.path.join(self._af3_dir, "command", self.lib_name)
        )

        self._af3_command_path = os.path.join(
            self._af3_command_subdir, f"run_af3_{self._gen_opt}.sh"
        )

        self._af3_struct_subdir = checkNgen_folder(
            os.path.join(self._af3_dir, "struct", self.lib_name)
        )

        self._full_abs_input_dir = os.path.abspath(self._af3_json_subdir)
        self._full_abs_output_dir = os.path.abspath(self._af3_struct_subdir)

        self._ifrerun = ifrerun

        self._sub_smiles = canonicalize_smiles(self.lib_info["substrate-smiles"])
        self._sub_dets = self.lib_info["substrate"]

        self._cofactor_smiles = canonicalize_smiles(
            ".".join(self.lib_info[f"{cofactor_dets}-smiles"])
        )
        self._cofactor_dets = "-".join(self.lib_info[cofactor_dets])

        self._joint_smiles = self._sub_smiles + "." + self._cofactor_smiles
        self._joint_dets = self._sub_dets + "_" + self._cofactor_dets

        self._full_model_path = full_model_path
        self._full_db_path = full_db_path

        self._gen_msa()

        self._make_var_json()
        self._write_var_exe()

    def _gen_msa(self):

        """
        A method to generate the MSA for each enzyme and then to be used for different variants.
        """

        # check if exist and ifrerun is False
        if os.path.exists(self._af3_msa_outpath) and not self._ifrerun:
            return

        json_data = {
            "name": f"{self.protein_name}",
            "sequences": [{"protein": {"id": "A", "sequence": f"{self.parent_seq}"}}],
            "modelSeeds": [1],
            "dialect": "alphafold3",
            "version": 1,
        }

        # Save the JSON to generate msa
        with open(self._af3_msa_inpath, "w") as json_file:
            json.dump(json_data, json_file, indent=4)

        # use the json for the docker run
        with open(self._af3_msa_exe, "w") as sh_file:
            sh_file.write("#!/bin/bash\n\n")

            docker_command = f"""docker run -it \\
                --volume {os.path.abspath(self._af3_msa_indir)}:/root/af_input \\
                --volume {os.path.abspath(self._af3_msa_dir)}:/root/af_output \\
                --volume {self._full_model_path}:/root/models \\
                --volume {self._full_db_path}:/root/public_databases \\
                --gpus all \\
                alphafold3 \\
                python run_alphafold.py \\
                --json_path=/root/af_input/{os.path.basename(self._af3_msa_inpath)} \\
                --jackhmmer_n_cpu=32 \\
                --nhmmer_n_cpu=32 \\
                --model_dir=/root/models \\
                --output_dir=/root/af_output \\
                --run_data_pipeline=True \\
                --run_inference=False\n\n"""
            sh_file.write(docker_command)

        # Make the shell script executable
        os.chmod(self._af3_msa_exe, 0o755)

        # Log file for the output
        log_file = os.path.splitext(self._af3_msa_exe)[0] + ".log"

        # Execute the generated shell script
        try:
            with open(log_file, "w") as log:
                subprocess.run(
                    ["sudo", self._af3_msa_exe],
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True,
                )
            print(f"MSA generation completed successfully. Logs saved to {log_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error during MSA generation: {e}")
            with open(log_file, "a") as log:
                log.write("\nDocker execution failed.\n")

    def _make_var_json(self):
        """
        A method to generate the AlphaFold3 JSON files for each variant.
        """

        msa_json = load_json(self._af3_msa_outpath)

        for (
            var,
            seq,
        ) in tqdm(self.df[[self._var_col_name, self._seq_col_name]].values):

            # Create the JSON structure
            json_data = {
                "name": f"{var}",
                "sequences": [
                    {
                        "protein": {
                            "id": "A",
                            "sequence": f"{seq}",
                            "modifications": msa_json.get(
                                "modifications", None
                            ),
                            "unpairedMsa": msa_json.get(
                                "unpairedMsa", None
                            ),
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

                for j, cofactor_dets, cofactor_smiles in enumerate(
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

    def _write_var_exe(self):

        """
        A method to write the shell script to execute the AlphaFold3 inference.
        """

        # TOOD: check if exist and ifrerun is False

        with open(self._af3_command_path, "w") as sh_file:
            sh_file.write("#!/bin/bash\n\n")

            # Iterate over all JSON files
            for json_file in os.listdir(self._af3_json_subdir):

                if json_file.endswith(".json"):

                    docker_command = f"""docker run -it \\
            --volume {self._full_abs_input_dir}:/root/af_input \\
            --volume {self._full_abs_output_dir}:/root/af_output \\
            --volume {self._full_model_path}:/root/models \\
            --volume {self._full_db_path}:/root/public_databases \\
            --gpus all \\
            alphafold3 \\
            python run_alphafold.py \\
            --json_path=/root/af_input/{json_file} \\
            --jackhmmer_n_cpu=32 \\
            --nhmmer_n_cpu=32 \\
            --output_dir=/root/af_output \\
            --run_data_pipeline=False \\
            --run_inference=True\n\n"""
                    sh_file.write(docker_command)

    @property
    def af3_dir(self) -> str:
        """Return the directory for the AlphaFold3 directory."""
        return self._af3_dir

    @property
    def af3_command_path(self) -> str:
        """
        Return the path for the shell script
        to execute the AlphaFold3 inference for the library.
        """
        return self._af3_command_path


def run_af3_prep(
    pattern: str | list = "data/meta/not_scaled/*.csv",
    gen_opt: str = "joint",
    kwargs: dict = {},
):
    """
    Run af3 for all libraries

    Args:
    - pattern: str | list: the pattern for the input csv files
    - gen_opt: str: The generation option for the AF3 structure
    - kwargs: dict: The arguments for the AF3Data class
    """

    all_exe_paths = []

    if isinstance(pattern, str):
        lib_list = sorted(glob(pattern))
    else:
        lib_list = deepcopy(pattern)

    for lib in lib_list:
        print(f"Running AF3 for {lib}...")
        AF3 = AF3Prep(input_csv=lib, gen_opt=gen_opt, **kwargs)
        all_exe_paths.append(AF3.af3_command_path)

    # get where this runner script will be saved
    runner_path = os.path.join(AF3.af3_dir, f"run_af3_{gen_opt}.sh")

    # write the runner script
    with open(runner_path, "w") as runner_file:
        runner_file.write("#!/bin/bash\n\n")

        for exe_path in all_exe_paths:
            relative_exe_path = exe_path.replace(
                os.path.normpath(AF3.af3_dir) + "/", ""
            )
            runner_file.write(f"sudo bash {relative_exe_path}\n")

    # Make the runner script executable
    os.chmod(runner_path, 0o755)

    # now use subprocess to run the runner script
    subprocess.run(["sudo", runner_path], check=True)