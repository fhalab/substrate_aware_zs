# File to automatically run LigandMPNN as a ZS
# Modify from Lukas Radtke {lradtke@caltech.edu}
# Date: 22/11/24

from __future__ import annotations

import subprocess
import re
import os
import math
from glob import glob
from tqdm import tqdm
from copy import deepcopy

import pandas as pd
import torch

from REVIVAL.preprocess import ZSData
from REVIVAL.util import checkNgen_folder


LigandMPNN_MODEL_LIST = deepcopy(
    [
        "ligandmpnn_v_32_005_25.pt",
        "ligandmpnn_v_32_010_25.pt",
        "ligandmpnn_v_32_020_25.pt",
        "ligandmpnn_v_32_030_25.pt",
    ]
)


LigandMPNN_MODEL_DICT = {int(m.split("_")[-2]): m for m in LigandMPNN_MODEL_LIST}


LigandMPNN_MODEL_PARAM_DIR = "/disk2/fli/LigandMPNN/model_params"

LigandMPNN_MAIN_DIR = "/disk2/fli/LigandMPNN"


class LigandmpnnData(ZSData):
    def __init__(
        self,
        input_csv: str,
        docked_struct_dir: str = "data/structure/docked",
        ligandmpnn_main_dir: str = LigandMPNN_MAIN_DIR,  # for running the .py file
        ligandmpnn_model_dir: str = LigandMPNN_MODEL_PARAM_DIR,  # for model .pt file
        noise_level: int = 20,  # [5, 10, 20, 30]
        scale_fit: str = "parent",
        chain_id: str = "A",  # the chain for enzyme to be modifed
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        mut_col_name: str = "mut",
        pos_col_name: str = "pos",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        output_dir: str = "zs/ligandmpnn",
        scoring_opt: str = "autoregressive",
        number_of_batches: int = 100,
        batch_size: int = 64,
    ):

        """
        Use AF3 or Chai docked structures

        The outputs inlucde:
            - The input pdb file under `input_struct` subfolder,
                converted from cif or made a copy
            - A .pt file saved under `raw_output` subfolder, with keys:
                ['logits', 'probs', 'log_probs', 'decoding_order',
                'native_sequence', 'mask', 'chain_mask', 'seed',
                'alphabet', 'residue_names', 'sequence', 'mean_of_probs', 'std_of_probs']
            - The corresponding .csv file extracting the ZS scores saved under the `scores` subfolder

        Args:
        - scoring_opt: str, default "autoregressive"
        """

        super().__init__(
            input_csv=input_csv,
            scale_fit=scale_fit,
            combo_col_name=combo_col_name,
            var_col_name=var_col_name,
            mut_col_name=mut_col_name,
            pos_col_name=mut_col_name,
            seq_col_name=seq_col_name,
            fit_col_name=fit_col_name,
        )

        self._output_dir = output_dir
        self._docked_struct_dir = docked_struct_dir
        self._chain_id = chain_id
        self._noise_level = noise_level
        self._scoring_opt = scoring_opt
        self._number_of_batches = number_of_batches
        self._batch_size = batch_size

        self._model_path = os.path.join(
            ligandmpnn_model_dir, LigandMPNN_MODEL_DICT[self._noise_level]
        )

        # get the score.py from the ligandmpnn path
        self._infernece_script_path = os.path.join(ligandmpnn_main_dir, "score.py")

        # preprocess cif to pdb
        self._preprocess_cif()

        print(f"Running model inference for {self.lib_name}")
        self._run_inference()

        print(f"Parsing score for {self.lib_name}")
        score_df = self._get_lib_score()

        # merge the score with the input csv fitness
        self._score_df = pd.merge(
            self.input_df[[self._combo_col_name, self._fit_col_name]],
            score_df,
            on=self._combo_col_name,
            how="left",
        )

        # save the score to csv
        self._score_df.to_csv(self.score_df_path, index=False)

    def _preprocess_cif(self) -> None:
        """
        Preprocess the cif file to pdb file
        """
        # check if the cif file exists
        if not os.path.exists(self.cif_path):
            raise FileNotFoundError(f"{self.cif_path} does not exist.")

        # check if the pdb file exists
        if not os.path.exists(self.pdb_path):
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

    def _run_inference(self) -> None:

        """
        Prepares the input for the model and runs ligandmpnn
        """

        cli_ligandmpnn = [
            "python",
            self._infernece_script_path,
            "--model_type",
            "ligand_mpnn",
            "--checkpoint_ligand_mpnn",
            self._model_path,
            "--seed",
            "12345",
            f"--{self._scoring_opt}_score",  # TODO write custom joint prob function
            "1",  # to obtain aa-level probabilities in the ouput file
            "--use_sequence",
            "1",  # autoregressively decodes residues to design, score depends on decoding order
            "--number_of_batches",
            str(self._number_of_batches),
            "--pdb_path",
            self.pdb_path,
            "--redesigned_residues",
            self.mutated_pos,  # takes the set of all mutated positions
            "--out_folder",
            self.raw_pt_dir,
            "--batch_size",
            str(self._batch_size),
        ]

        try:
            subprocess.run(cli_ligandmpnn, capture_output=True, text=True)
        except Exception as e:
            print(f"Error running ligandmpnn: {e}")

    def _get_lib_score(self) -> list:

        """
        Extract the ZS scores from the output file
        """

        score_df = []

        # load .pt file
        raw_pt = torch.load(self.raw_pt_path)

        for combo in self.input_df[self._combo_col_name].to_list():

            log_likelyhoods = 0

            for mutated_position, combo_aa in zip(self.mutated_pos_list, combo):
                log_likelyhoods += math.log(
                    raw_pt["mean_of_probs"][mutated_position][combo_aa]
                )

            score_df.append(
                {
                    self._combo_col_name: combo,
                    "ligandmpnn_score": log_likelyhoods,
                }
            )

        
        return pd.DataFrame(score_df)

    
    @property
    def mutated_pos_list(self) -> list:
        """
        Return a list of mutations for ligandmpnn: [A165, A183, A301]
        """
        return [str(self._chain_id) + str(pos) for pos in self.mut_pos_list]
    
    @property
    def mutated_pos(self) -> str:
        """
        Return a set of mutations for ligandmpnn: 'A165 A183 A301'
        """
        return " ".join(self.mutated_pos_list)

    @property
    def ligandmpnn_dets(self) -> str:
        """
        Return the details of the ligandmpnn run
        """
        return f"{self._scoring_opt}_{self._noise_level}"

    @property
    def raw_pt_dir(self) -> str:
        """
        Path for the ligandmpnn .pt file output under the `raw_output` subfolder.
        """
        return checkNgen_folder(
            os.path.join(self._output_dir, "raw_output", self.ligandmpnn_dets)
        )

    @property
    def raw_pt_path(self) -> str:
        """
        Filepath to the corresponding .pt file for a given .csv of same name.
        """
        return os.path.join(self.raw_pt_dir, f"{self.lib_name}.pt")

    @property
    def cif_path(self) -> str:
        """
        Filepath to the corresponding .cif file for a given .csv of same name.
        """
        return os.path.join(self._docked_struct_dir, f"{self.lib_name}.cif")

    @property
    def pdb_path(self) -> str:
        """
        Filepath to the corresponding .pdb file for a given .csv of same name.
        """
        input_pdb_dir = checkNgen_folder(os.path.join(self._output_dir, "input_struct"))

        return os.path.join(input_pdb_dir, f"{self.lib_name}.pdb")

    @property
    def score_dir(self) -> str:
        """
        Path for the ligandmpnn .csv file output under the `scores` subfolder.
        """
        return checkNgen_folder(
            os.path.join(self._output_dir, "scores", self.ligandmpnn_dets)
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


def run_ligandmpnn(pattern: str | list = None, kwargs: dict = {}):

    if isinstance(pattern, str):
        lib_list = sorted(glob(pattern))
    else:
        lib_list = deepcopy(pattern)

    for lib in tqdm(lib_list):
        print(f"Running LigandMPNN for {lib}...")
        LigandmpnnData(input_csv=lib, **kwargs)