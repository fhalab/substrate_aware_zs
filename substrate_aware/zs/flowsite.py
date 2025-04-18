# File to automatically run Flowsite as a ZS
# Modified from Lukas Radtke {lradtke@caltech.edu}
# Date: 06/12/24

import os

import subprocess
from glob import glob
from tqdm import tqdm
from copy import deepcopy

import numpy as np
import pandas as pd


from substrate_aware.preprocess import ZSData
from substrate_aware.global_param import ENZYME_INFO_DICT
from substrate_aware.util import checkNgen_folder, calculate_ligand_centroid


FLOWSITE_MAIN_DIR = "/disk2/fli/FlowSite"
FLOWSITE_MODEL_PARAM_DIR = "/disk2/fli/FlowSite/model_params"

# residues_canonical_1letter in https://github.com/HannesStark/FlowSite/blob/main/utils/featurize.py
# omit the alst one misc
FLOWSITE_AA_ORDER = [
    "A",
    "R",
    "N",
    "D",
    "C",
    "Q",
    "E",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V",
]
FLOWSITE_AA_IDX = {aa: i for i, aa in enumerate(FLOWSITE_AA_ORDER)}

# model and hyperparameter choices according to https://github.com/HannesStark/FlowSite
FLOWSITE_MODEL_DETS = {
    1: {
        "checkpoint": "pocket_gen/lf5t55w4/checkpoints/best.ckpt",
        "ns": 48,
        "nv": 12,
    },  # from csv inference example
    2: {
        "checkpoint": "pocket_gen/b1ribx1a/checkpoints/best.ckpt",
        "ns": 32,
        "nv": 8,
    },  # recommended for pocket_def_residues or pocket_def_center
}

FLOWSITE_DEFAULT_PARAMS = {
    "num_inference": 100,
    "batch_size": 64,
    "fake_constant_dur": 100000,
    "fake_decay_dur": 10,
    "fake_ratio_start": 0.2,
    "fake_ratio_end": 0.2,
    "residue_loss_weight": 0.2,
    "flow_matching_sigma": 0.5,
    "prior_scale": 1,
    "num_angle_pred": 5,
    "pocket_residue_cutoff": 12,
}


class FlowsiteData(ZSData):
    def __init__(
        self,
        input_csv=str,
        combo_col_name: str = "AAs",
        fit_col_name: str = "fitness",
        chain_id: str = "A",
        dock_opt: str = "substrate",  # 'substrate' or 'joint'
        cofactor_dets: str = "cofactor",  # 'cofactor' or 'inactivated-cofactor'
        in_struct_dir: str = "data/structure/docked",
        output_dir: str = "zs/flowsite",
        flowsite_inference_opt: str = "pocket_def_residues",  # or pocket_def_center
        flowsite_model_opt: int = 2,
        flowsite_main_dir: str = FLOWSITE_MAIN_DIR,  # for running the .py file
        flowsite_model_dir: str = FLOWSITE_MODEL_PARAM_DIR,  # for checkpoint .ckpt file
        flowsite_model_dets: dict = FLOWSITE_MODEL_DETS,
        flowsite_params: dict = FLOWSITE_DEFAULT_PARAMS,
        regen: bool = False,
    ):

        """
        Use closest pdb, AF3 or Chai docked structures

        The outputs inlucde:
            - The input pdb file under `input_struct` subfolder,
                converted from cif or made a copy
            - A .npy file saved under `raw_output/lib_name/complexid0` subfolder, named `desgined_logits`
                with deminsions: (num_inference, num_residues, 21)
            - The corresponding .csv file extracting the ZS scores saved under the `scores` subfolder
        """

        super().__init__(
            input_csv=input_csv,
            combo_col_name=combo_col_name,
            fit_col_name=fit_col_name,
        )

        self._chain_id = chain_id
        self._dock_opt = dock_opt
        self._cofactor_dets = cofactor_dets

        self._in_struct_dir = in_struct_dir
        self._output_dir = output_dir

        self._inference_path = os.path.join(flowsite_main_dir, "inference.py")

        self._model_dets = flowsite_model_dets[flowsite_model_opt]
        self._checkpoint_path = os.path.join(
            flowsite_model_dir, self._model_dets["checkpoint"]
        )

        self._flowsite_params = flowsite_params
        self._flowsite_inference_opt = flowsite_inference_opt

        self._flowsite_dets = f"{self._flowsite_inference_opt}_model{flowsite_model_opt}-{dock_opt}-{cofactor_dets}"

        if not regen and os.path.exists(self.raw_npy_path):
            print(
                f"Flowsite npy output already exists for {self.lib_name} and not regen."
            )

        else:
            self._preprocess_struct()

            print(f"Running flowsite for {self.lib_name}...")
            self._run_flowsite()

        print(f"Creating score .csv for {self.lib_name}...")
        score_df = self._get_flowsite_score()

        # merge the score_df with the input_df
        self._score_df = pd.merge(
            self.df[[self._combo_col_name, self._fit_col_name]],
            score_df,
            on=self._combo_col_name,
            how="left",
        )

        # save the score_df to a csv file
        self._score_df.to_csv(self.score_df_path, index=False)

    def _preprocess_struct(self) -> None:
        """
        Preprocess the cif file to pdb file or make a copy of the pdb file
        """

        # check if the pdb file exists
        if not os.path.exists(self.pdb_path):

            if os.path.exists(self.cif_path.replace(".cif", ".pdb")):
                # copy the pdb file to the input_struct folder
                subprocess.run(
                    [
                        "cp",
                        self.cif_path.replace(".cif", ".pdb"),
                        self.pdb_path,
                    ],
                    capture_output=True,
                )
            else:
                # convert the cif to pdb
                subprocess.run(
                    [
                        "obabel",
                        self.cif_path,
                        "-O",
                        self.pdb_path,
                    ],
                    capture_output=True,
                )

    def _run_flowsite(self) -> None:

        """
        Prepares the inputs for the model and runs FlowSite
        """

        # create a smiles for substrate and possibly cofactor
        if self._dock_opt == "substrate":
            smiles = self.substrate_smiles
        elif self._dock_opt == "joint":
            smiles = self.joint_smiles
        else:
            raise ValueError(f"Invalid docking option: {self.dock_opt}")

        # define console command and execute subprocess
        cli_cmd = [
            "python",
            os.path.abspath(self._inference_path),
            "--num_inference",
            str(self._flowsite_params["num_inference"]),
            "--batch_size",
            str(self._flowsite_params["batch_size"]),
            "--out_dir",
            os.path.abspath(self.raw_npy_dir),
            "--protein",
            os.path.abspath(self.pdb_path),
            "--smiles",
            smiles,
            "--design_residues",
            str(self.mutated_pos),
            "--checkpoint",
            os.path.abspath(self._checkpoint_path),
            "--run_test",
            "--run_name",
            "inference1",
            "--layer_norm",
            "--fake_constant_dur",
            str(self._flowsite_params["fake_constant_dur"]),
            "--fake_decay_dur",
            str(self._flowsite_params["fake_decay_dur"]),
            "--fake_ratio_start",
            str(self._flowsite_params["fake_ratio_start"]),
            "--fake_ratio_end",
            str(self._flowsite_params["fake_ratio_end"]),
            "--residue_loss_weight",
            str(self._flowsite_params["residue_loss_weight"]),
            "--use_tfn",
            "--use_inv",
            "--time_condition_tfn",
            "--correct_time_condition",
            "--time_condition_inv",
            "--time_condition_repeat",
            "--flow_matching",
            "--flow_matching_sigma",
            str(self._flowsite_params["flow_matching_sigma"]),
            "--prior_scale",
            str(self._flowsite_params["prior_scale"]),
            "--num_angle_pred",
            str(self._flowsite_params["num_angle_pred"]),
            "--ns",
            str(self._model_dets["ns"]),
            "--nv",
            str(self._model_dets["nv"]),
            "--batch_norm",
            "--self_condition_x",
            "--self_condition_inv",
            "--no_tfn_self_condition_inv",
            "--self_fancy_init",
            "--pocket_residue_cutoff",
            str(self._flowsite_params["pocket_residue_cutoff"]),
        ]

        if self._flowsite_inference_opt == "pocket_def_center":
            cli_cmd.extend(
                [
                    "--pocket_def_center",
                    self.ligand_centroid,
                ]
            )
        else:
            cli_cmd.extend(
                [
                    "--pocket_def_residues",
                    str(self.mutated_pos),
                ]
            )

        process = subprocess.run(cli_cmd, capture_output=True, text=True)

        # Print the outputs
        print(
            "Return Code:", process.returncode
        )  # Check if the process exited successfully
        print("Standard Output:", process.stdout)  # Output of the command
        print("Standard Error:", process.stderr)  # Error messages, if any

        # Handle errors explicitly
        if process.returncode != 0:
            print("An error occurred while executing the command.")
            print("Error details:", process.stderr)

    def _get_flowsite_score(self) -> pd.DataFrame:

        """
        Extracts the ZS scores from flowsite output files and writes them to a csv.
        """

        score_df = []

        for combo in self.df[self._combo_col_name]:
            combo_scores = []
            for rep in range(self._flowsite_params["num_inference"]):
                raw_score = self.designed_logits[rep]
                combo_scores.append(
                    sum(
                        [
                            raw_score[i][FLOWSITE_AA_IDX[res]]
                            for i, res in enumerate(combo)
                        ]
                    )
                )
            score_df.append(
                {
                    self._combo_col_name: combo,
                    "flowsite_score": np.mean(combo_scores),
                    "flowsite_std": np.std(combo_scores),
                }
            )

        return pd.DataFrame(score_df)

    @property
    def mutated_pos(self) -> str:
        """
        Extracts the mutations position and prepares them in the expected format of flowsite (e.g. 'A165, A183, A301').
        Assumes the same mutation positions within the entire campaign.
        """
        return ",".join([self._chain_id + str(pos) for pos in self.mut_pos_list])

    @property
    def substrate_smiles(self) -> str:
        """
        Extracts the substrate smiles from the input csv.
        Assumes the same substrate within the entire campaign.
        """
        return self.lib_info["substrate-smiles"]

    @property
    def joint_smiles(self) -> str:
        """
        Extracts the substrate and cofactor smiles from the input csv.
        Assumes the same substrate and cofactor within the entire campaign.
        """
        return self.lib_info["substrate-smiles"] + ".".join(
            self.lib_info[f"{self._cofactor_dets}-smiles"]
        )

    @property
    def ligand_centroid(self) -> np.ndarray:
        """
        Extracts the ligand centroid from the input pdb file.
        """
        return ",".join(map(str,calculate_ligand_centroid(
            pdb_file=self.pdb_path,
            ligand_info=ENZYME_INFO_DICT[self.protein_name]["ligand-info"]
        )))

    @property
    def raw_npy_dir(self) -> str:
        """
        Path for the output .npy file under the `raw_output` subfolder.
        """
        return checkNgen_folder(
            os.path.join(
                self._output_dir, "raw_output", self._flowsite_dets, self.lib_name
            )
        )

    @property
    def raw_npy_path(self) -> str:
        """
        Filepath to the corresponding .npy file for a given .csv of same name.
        """
        return os.path.join(self.raw_npy_dir, "complexid0/designed_logits.npy")

    @property
    def designed_logits(self) -> np.ndarray:
        """
        Return the logits array
        """
        return np.load(self.raw_npy_path)

    @property
    def cif_path(self) -> str:
        """
        Filepath to the corresponding .cif file for a given .csv of same name.
        """
        return os.path.join(self._in_struct_dir, f"{self.lib_name}.cif")

    @property
    def pdb_path(self) -> str:
        """
        Filepath to the corresponding .pdb file for a given .csv of same name.
        """
        input_pdb_dir = checkNgen_folder(
            os.path.join(self._output_dir, "input_struct", self._flowsite_dets)
        )

        return os.path.join(input_pdb_dir, f"{self.lib_name}.pdb")

    @property
    def score_dir(self) -> str:
        """
        Path for the flowsite .csv file output under the `scores` subfolder.
        """
        return checkNgen_folder(
            os.path.join(self._output_dir, "scores", self._flowsite_dets)
        )

    @property
    def score_df_path(self) -> str:
        """
        Filepath to the corresponding .csv file for a given .csv of same name.
        """
        return os.path.join(self.score_dir, f"{self.lib_name}.csv")

    @property
    def score_df(self) -> pd.DataFrame:
        """
        Return the score dataframe
        """
        return self._score_df


def run_flowsite(
    pattern: str | list = None,
    flowsite_inference_opt: str = "pocket_def_residues",
    flowsite_model_opt: int = 2,
    kwargs: dict = {}):

    if isinstance(pattern, str):
        lib_list = sorted(glob(pattern))
    else:
        lib_list = deepcopy(pattern)

    for lib in tqdm(lib_list):
        FlowsiteData(
            input_csv=lib,
            flowsite_inference_opt=flowsite_inference_opt,
            flowsite_model_opt=flowsite_model_opt,
            **kwargs
        )