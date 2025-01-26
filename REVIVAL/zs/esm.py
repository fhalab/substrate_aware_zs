"""
A script for generating EMS zs scores
"""

from __future__ import annotations

import os
from glob import glob
from copy import deepcopy
from tqdm import tqdm

import numpy as np
import pandas as pd

import esm
import torch

from REVIVAL.preprocess import ZSData
from REVIVAL.util import checkNgen_folder


class ESM(ZSData):

    """
    A class for generating ESM scores
    """

    def __init__(
        self,
        input_csv: str,
        esm_model_name: str = "esm2_t33_650M_UR50D",
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        mut_col_name: str = "mut",
        pos_col_name: str = "pos",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        seq_dir: str = "data/seq",
        zs_dir: str = "zs",
        esm_dir: str = "esm",
        regen_esm=False,
    ):

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

        self._esm_model_name = esm_model_name

        (
            self._model,
            self._alphabet,
            self._batch_converter,
            self._device,
        ) = self._infer_model()

        self._mask_string, self._cls_string, self._eos_string = (
            self._alphabet.mask_idx,
            self._alphabet.cls_idx,
            self._alphabet.eos_idx,
        )

        self._alphabet_size = len(self._alphabet)

        self._esm_dir = checkNgen_folder(os.path.join(zs_dir, esm_dir))
        self._esm_output_dir = checkNgen_folder(os.path.join(self._esm_dir, "output"))
        self._esm_path = os.path.join(self._esm_output_dir, f"{self.lib_name}.csv")

        self._logit_dir = checkNgen_folder(
            os.path.join(self._esm_dir, "logits", self._esm_model_name)
        )
        self._logit_path = os.path.join(self._logit_dir, f"{self.lib_name}.npy")

        if (
            self._logit_path != ""
            and os.path.exists(self._logit_path)
            and not (regen_esm)
        ):
            print(f"{self._logit_path} exists and regen_esm = {regen_esm}. Loading...")
            self._logits = np.load(self._logit_path)
        else:
            print(f"Generating {self._logit_path}...")
            self._logits = self._get_logits()

        # generate esm score
        self._esm_df = self._run_esm()

        # if the esm score for the WT is empty then update the esm score to 0

        # save the esm score df
        self._esm_df.to_csv(self._esm_path, index=False)

    def _infer_model(self):

        """
        Infer the model
        """

        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        esm_pretrained = getattr(esm.pretrained, self._esm_model_name)
        model, alphabet = esm_pretrained()
        batch_converter = alphabet.get_batch_converter()
        model.eval()
        model = model.to(device)
        print("Using device:", device)
        return model, alphabet, batch_converter, device

    def _get_logits(self):

        """
        Get the logits for the wild type sequence
        """

        data_wt = [("WT", self.parent_seq)]
        # Get Batch tokens for WT data
        batch_labels_wt, batch_strs_wt, batch_tokens_wt = self._batch_converter(data_wt)

        logits = np.zeros((len(self.parent_seq), self._alphabet_size))

        for (i, seq) in enumerate(data_wt[0][1]):
            batch_tokens_masked = batch_tokens_wt.clone()
            batch_tokens_masked[0, i] = self._alphabet.mask_idx
            batch_tokens_masked = batch_tokens_masked.to(self._device)

            with torch.no_grad():
                token_probs = (
                    torch.log_softmax(
                        self._model(batch_tokens_masked)["logits"], dim=-1
                    )
                    .cpu()
                    .numpy()
                )

            logits[i] = token_probs[0, i + 1]

        # save the logits
        np.save(self._logit_path, logits)

        return logits

    def _get_mutant_prob(self, mt, wt, pos):

        """
        Get the probability of the mutant given the wild type sequence at certain position.
        """

        wt_idx = self._alphabet.get_idx(wt)
        mt_idx = self._alphabet.get_idx(mt)

        return self._logits[pos, mt_idx] - self._logits[pos, wt_idx]

    def _get_esm_score(self, df):

        """
        Run ESM model for all variants in the data set

        Input:  - logits: Logits of the wild type sequence
                - df: Data set containing the variants, loops trough column = 'mut' and 'pos'
                - _sum: If True, the sum of the probabilities is calculated.
                    If False, the mean of the probabilities is calculated

        Output: - Score for each variant
        """

        score = np.zeros(len(df))
        wt_sequence = list(self.parent_seq)

        for i, (mut_list, pos_list) in enumerate(
            zip(df[self._mut_col_name], df[self._pos_col_name])
        ):
            s = np.zeros(len(mut_list))
            for j, mt in enumerate(mut_list):
                if mt == "WT":
                    score[i] = 0
                    continue

                elif mt == "NA":
                    score[i] = np.nan
                    continue

                else:
                    # Position of the mutation with python indexing
                    pos = int(pos_list[j]) - 1
                    wt = wt_sequence[pos]
                    s[j] = self._get_mutant_prob(mt=mt, wt=wt, pos=pos)

                score[i] += s.sum()

        return score

    def _run_esm(self):

        """
        Get any score for each variant in the data set
        """

        df_n_list = []

        # TODO test
        wt_df = self.df[self.df[self._var_col_name] == "WT"].copy()
        # add if wt_df is not empty
        if not wt_df.empty:
            wt_df.loc[:, "esm_score"] = 0

        df_n_list.append(wt_df)

        # Get the n mutant scores
        for i in tqdm(list(range(self.max_n_mut + 1))[1:]):
            # Filter out n mutants
            df_n = self.df[self.df["n_mut"] == i].copy()
            if df_n.empty:  # Check if the DataFrame is empty after filtering
                assert "Data set is empty"
                continue

            score_n = self._get_esm_score(df_n)

            # Add column with number of mutations
            df_n.loc[:, "esm_score"] = score_n

            df_n_list.append(df_n)

        return pd.concat(df_n_list, axis=0)


def run_all_esm(
    pattern: str | list = "data/meta/not_scaled/*",
    esm_model_name: str = "esm2_t33_650M_UR50D",
    combo_col_name: str = "AAs",
    var_col_name: str = "var",
    mut_col_name: str = "mut",
    pos_col_name: str = "pos",
    seq_col_name: str = "seq",
    fit_col_name: str = "fitness",
    seq_dir: str = "data/seq",
    zs_dir: str = "zs",
    esm_dir: str = "esm",
    regen_esm=False,
):

    """
    Run all ESM scores
    """

    if isinstance(pattern, str):
        path_list = glob(pattern)
    else:
        path_list = deepcopy(pattern)

    for p in tqdm(path_list):

        print(f"Running ESM for {p}...")

        ESM(
            input_csv=p,
            esm_model_name=esm_model_name,
            combo_col_name=combo_col_name,
            var_col_name=var_col_name,
            mut_col_name=mut_col_name,
            pos_col_name=pos_col_name,
            seq_col_name=seq_col_name,
            fit_col_name=fit_col_name,
            seq_dir=seq_dir,
            zs_dir=zs_dir,
            esm_dir=esm_dir,
            regen_esm=regen_esm,
        )