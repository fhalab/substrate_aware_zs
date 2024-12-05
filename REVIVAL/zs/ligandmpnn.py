# TODO: modify and test the following

# File to automatically run LigandMPNN as a ZS
# Author: Lukas Radtke {lradtke@caltech.edu}
# Date: 22/11/24

import subprocess
import os
import math
from glob import glob
from copy import deepcopy

import pandas as pd
import re
import torch
from tqdm import tqdm


from REVIVAL.preprocess import ZSData
from REVIVAL.util import checkNgen_folder


class LigandmpnnData(ZSData):
    def __init__(
        self,
        input_csv: str,
        scale_fit: str = "parent",
        chain_id: str = "A",
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        mut_col_name: str = "mut",
        pos_col_name: str = "pos",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        structure_dir: str = "data/structure",
        ligandmpnn_dir: str = "../LigandMPNN",
        results_dir: str = "results",
        noise_levels: list = [5, 10, 20, 30],
    ):

        super().__init__(
            input_csv=input_csv,
            scale_fit=scale_fit,
            combo_col_name=combo_col_name,
            var_col_name=var_col_name,
            mut_col_name=mut_col_name,
            pos_col_name=mut_col_name,
            seq_col_name=seq_col_name,
            fit_col_name=fit_col_name,
            structure_dir=structure_dir,
        )

        self.ligandmpnn_dir = ligandmpnn_dir
        self.results_dir = results_dir
        self.chain_id = chain_id
        self.noise_levels = noise_levels
        self.structure_dir = structure_dir

        self.noise_dict = {
            5: "ligandmpnn_v_32_005_25.pt",
            10: "ligandmpnn_v_32_010_25.pt",
            20: "ligandmpnn_v_32_020_25.pt",
            30: "ligandmpnn_v_32_030_25.pt",
        }

        self.ligandmpnn_noise_models = [self.noise_dict[m] for m in self.noise_levels]

        self.inference_path = os.path.join(self.ligandmpnn_dir, "score.py")

        print(f"Running model inference for {self.lib_name}")
        self.inference()

        print(f"Creating score .csv for {self.lib_name}")
        self.score_mut_file()

    def extract_mutated_pos(self, csv) -> str:
        """
        Takes all var strings e.g. I165A:I183A:Y301V from a csv
        and converts them to a set of mutations for ligandmpnn: 'A165 A183 A301'
        """
        all_mutation = []
        for i in range(0, csv.shape[0]):
            var = csv.iloc[i][self._var_col_name]
            mutated_pos = re.findall(r"\d+", var)
            mutated_pos = [str(self.chain_id) + str(num) for num in mutated_pos]
            all_mutation.append(mutated_pos)
        all_mutation = sorted(set([mut for var in all_mutation for mut in var]))
        var_ligandmpnn_format = " ".join(all_mutation)
        return var_ligandmpnn_format

    def extract_noise_level_num(self, noise_level) -> str:
        """
        Converts noise-level model path to number
        """
        return noise_level[-9:-6]

    def inference(self) -> None:
        """
        Prepares the input for the model and runs ligandmpnn
        """
        csv = pd.read_csv(self._input_csv, keep_default_na=False)
        csv = csv[
            ~csv.apply(lambda row: row.astype(str).str.contains(r"\*").any(), axis=1)
        ].reset_index(
            drop=True
        )  # filters for stop codons, marked as '*'
        redesigned_residues = self.extract_mutated_pos(csv)

        for noise_level in self.ligandmpnn_noise_models:

            model_path = os.path.join(self.ligandmpnn_dir, "model_params", noise_level)
            self.model_results_path = os.path.join(
                self.results_dir,
                f"{os.path.basename(self._input_csv)[:-4]}_noise_level_{self.extract_noise_level_num(noise_level)}",
            )

            cli_ligandmpnn = [
                "python",
                self.inference_path,
                "--model_type",
                "ligand_mpnn",
                "--checkpoint_ligand_mpnn",
                model_path,
                "--seed",
                "12345",
                "--autoregressive_score",
                "1",  # to obtain aa-level probabilities in the ouput file
                "--use_sequence",
                "1",  # autoregressively decodes residues to design, score depends on decoding order, so I set n_batches high
                "--number_of_batches",
                "100",
                "--pdb_path",
                self.pdb_path,
                "--redesigned_residues",
                redesigned_residues,  # takes the set of all mutated positions in a csv and redesignes all of them for each variant, for ZS only the mutated positions of the variant are considered
                "--out_folder",
                self.model_results_path,
                "--batch_size",
                "1",
            ]

            print(f"Model noise_level: {noise_level[-9:-6]}")
            cli_return = subprocess.run(cli_ligandmpnn, capture_output=True, text=True)
            # print(cli_return)

    def score_mut_file(self) -> None:
        """
        A function for extracting the ZS scores from LigandMPNN output files and writing them to a csv.
        """
        score_df = pd.DataFrame(columns=["variant", "gt_fitness", "zs_score"])
        csv = pd.read_csv(self._input_csv, keep_default_na=False)
        csv = csv[
            ~csv.apply(lambda row: row.astype(str).str.contains(r"\*").any(), axis=1)
        ].reset_index(
            drop=True
        )  # filters for stop codons, marked as '*'

        for noise_level in self.ligandmpnn_noise_models:

            noise_level_num = self.extract_noise_level_num(noise_level)
            filename = os.path.basename(self._input_csv).replace(".csv", ".pt")
            filepath = os.path.join(
                self.results_dir,
                f"{filename[:-3]}_noise_level_{noise_level_num}",
                filename,
            )

            ZS = torch.load(filepath)

            rows = []
            for _, variant in tqdm(
                csv.iterrows(),
                desc=f"Extracing ZS scores for {self._input_csv} and noise {self.extract_noise_level_num(noise_level)}",
            ):

                var_pos = re.findall(r"\d+", variant[self._var_col_name])

                log_likelyhoods = []
                for i, mutated_position in enumerate(var_pos):
                    log_score = math.log(
                        ZS["mean_of_probs"][f"A{mutated_position}"][
                            variant[self._combo_col_name][i]
                        ]
                    )
                    log_likelyhoods.append(log_score)
                log_zs_score = sum(log_likelyhoods)
                rows.append(
                    {
                        "variant": variant[self._var_col_name],
                        "gt_fitness": variant[self._fit_col_name],
                        "zs_score": log_zs_score,
                    }
                )
                score_df = pd.DataFrame(rows)

            # saving comparison df as csv
            score_df_path = os.path.join(
                self.results_dir,
                "scores",
                f"{os.path.basename(self._input_csv)[:-4]}_noise_level_{noise_level_num}.csv",
            )
            checkNgen_folder(score_df_path)
            score_df.to_csv(score_df_path)
            print(f"ZS score .csv file saved to {score_df_path}.")

    @property
    def pdb_path(self) -> str:
        """
        Filepath to the corresponding .pdb file for a given .csv of same name.
        """
        return os.path.join(
            self._structure_dir,
            os.path.basename(self._input_csv).replace(".csv", ".pdb"),
        )


def run_ligandmpnn(pattern: str | list = None, kwargs: dict = {}):

    if isinstance(pattern, str):
        lib_list = glob(pattern)
    else:
        lib_list = deepcopy(pattern)

    for lib in lib_list:
        print(f"Running LigandMPNN for {lib}...")
        LigandmpnnData(input_csv=lib, **kwargs)