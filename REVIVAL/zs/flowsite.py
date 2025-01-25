# File to automatically run Flowsite as a ZS
# Modified from Lukas Radtke {lradtke@caltech.edu}
# Date: 06/12/24

import os
import re

import subprocess
from glob import glob
# from tqdm import tqdm
from copy import deepcopy

import numpy as np
import pandas as pd

from rdkit import Chem

# import torch

from REVIVAL.preprocess import ZSData
from REVIVAL.util import checkNgen_folder, add_hydrogens_to_smiles


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
    "num_inference": 10,
    "batch_size": 10,
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
        scale_fit: str = "parent",
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
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
    ):

        """
        Use closest pdb, AF3 or Chai docked structures

        The outputs inlucde:
            - The input pdb file under `input_struct` subfolder,
                converted from cif or made a copy
            - A .npy file saved under `raw_output` subfolder, with deminsions:
                (1, num_residues, 21)
            - The corresponding .csv file extracting the ZS scores saved under the `scores` subfolder
        """

        super().__init__(
            input_csv=input_csv,
            scale_fit=scale_fit,
            combo_col_name=combo_col_name,
            var_col_name=var_col_name,
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

        self._flowsite_dets = f"{flowsite_inference_opt}_model{flowsite_model_opt}-{dock_opt}-{cofactor_dets}"

        self._preprocess_struct()

        print(f"Running flowsite for {self.lib_name}...")
        self._run_flowsite()

        print(f"Creating score .csv for {self.lib_name}...")
        # self.score_mut_file()

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

        # try: 
        #     smiles = add_hydrogens_to_smiles(smiles)
        # except:
        #     pass

        # set up run device
        # cuda_available = torch.cuda.is_available()
        # env = os.environ.copy()
        # if cuda_available:
        #     env['CUDA_VISIBLE_DEVICES'] = '0'
        # else:
        #     env.pop('CUDA_VISIBLE_DEVICES', None)

        # define console command and execute subprocess
        cli_cmd = [
            "python",
            os.path.abspath(self._inference_path),
            "--num_inference",
            str(self._flowsite_params["num_inference"]),
            "--batch_size",
            str(self._flowsite_params["batch_size"]),
            "--out_dir",
            os.path.abspath(self.raw_npz_dir),
            "--protein",
            os.path.abspath(self.pdb_path),
            "--pocket_def_residues",
            str(self.mutated_pos),
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
        print(cli_cmd)
        process = subprocess.run(cli_cmd, capture_output=True, text=True)  # , env=env)

        # Print the outputs
        print("Return Code:", process.returncode)  # Check if the process exited successfully
        print("Standard Output:", process.stdout)  # Output of the command
        print("Standard Error:", process.stderr)  # Error messages, if any

        # Handle errors explicitly
        if process.returncode != 0:
            print("An error occurred while executing the command.")
            print("Error details:", process.stderr)

    # def score_mut_file(self) -> None:
    #     """
    #     A function for extracting the ZS scores from flowsite output files and writing them to a csv.
    #     """

    #     score_df = pd.DataFrame(columns=["variant", "gt_fitness", "zs_score"])
    #     csv = pd.read_csv(self._input_csv, keep_default_na=False)
    #     csv = csv[~csv.apply(lambda row: row.astype(str).str.contains(r"\*").any(), axis=1)].reset_index(drop=True)  # filters for stop codons, marked as '*'

    #     enzyme_name = csv['enzyme'][1] # assumes its the same enzyme for every variant
    #     substrate_name = csv['substrate'][1] # assumes its the same enzyme for every variant
    #     output_dir = os.path.join(self.results_dir, f'{enzyme_name}-{substrate_name}_inference', 'complexid0', 'designed_logits.npy')

    #     # load flowsite output
    #     ZS = np.load(output_dir)

    #     # extract ZS scores from flowsite output
    #     rows = []
    #     for _, variant in tqdm(csv.iterrows(), desc=f'Extracing ZS scores for {self._input_csv}'):
    #         log_likelyhoods_avg = []
    #         for i, res in enumerate(variant[self._combo_col_name]):
    #             log_score = sum([ZS[inference][i][FLOWSITE_AA_IDX[res]] for inference in range(self.num_inference)])
    #             log_likelyhoods_avg.append(log_score)
    #         log_zs_score = sum(log_likelyhoods_avg)
    #         rows.append({'variant': variant[self._var_col_name],
    #                      'gt_fitness': variant[self._fit_col_name],
    #                      'zs_score': log_zs_score})
    #     score_df = pd.DataFrame(rows)
    #     print(score_df['gt_fitness'].corr(score_df['zs_score'], method='spearman'))

    #     # save score_df
    #     score_df_path = os.path.join(self.results_dir, f'{os.path.basename(self._input_csv)[:-4]}_score_df.csv')
    #     score_df.to_csv(score_df_path)
    #     print(f'ZS score .csv file saved to {score_df_path}.')

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
    def raw_npz_dir(self) -> str:
        """
        Path for the output .npy file under the `raw_output` subfolder.
        """
        return checkNgen_folder(
            os.path.join(self._output_dir, "raw_output", self._flowsite_dets)
        )

    @property
    def raw_npz_path(self) -> str:
        """
        Filepath to the corresponding .pt file for a given .csv of same name.
        """
        return os.path.join(self.raw_npz_dir, f"{self.lib_name}.npz")

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
        input_pdb_dir = checkNgen_folder(os.path.join(self._output_dir, "input_struct", self._flowsite_dets))

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


def run_flowsite(pattern: str | list = None, kwargs: dict = {}):

    if isinstance(pattern, str):
        lib_list = glob(pattern)
    else:
        lib_list = deepcopy(pattern)

    for lib in lib_list:
        FlowsiteData(input_csv=lib, **kwargs)