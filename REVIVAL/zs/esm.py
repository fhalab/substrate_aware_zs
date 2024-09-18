"""
A script for generating EMS zs scores
"""

from __future__ import annotations

# Import packages
import os
import numpy as np
import pandas as pd

# ESM
import esm

# from esm.model.msa_transformer import MSATransformer

# Pytorch
import torch

from REVIVAL.global_param import LIB_INFO_DICT
from REVIVAL.preprocess import LibData
from REVIVAL.util import checkNgen_folder


class ZSData(LibData):
    """
    A class for generating ZS scores
    """
    def __init__(self,
        input_csv: str,
        scale_fit: str,
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        mut_col_name: str = "mut",
        pos_col_name: str = "pos",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        seq_dir: str = "data/seq"
        ):

        """
        - mut_col_name, str: the column name for the mutations
            ie ['A', 'D']
        - pos_col_name, str: the column name for the positions
            ie [39, 40]
        """

        super().__init__(input_csv, scale_fit, combo_col_name, var_col_name, seq_col_name, fit_col_name, seq_dir)

        self._mut_col_name = mut_col_name
        self._pos_col_name = pos_col_name


    def _append_mut_dets(self, combo: str) -> str:

        """
        Append mut details from the combo column

        Args:
        - combo, str: the variants sequence

        Returns:
        - list: the list of mutated AA
        - list: the list of mutated positions
        """

        mut_list = []
        pos_list = []

        for i, (mut, wt) in enumerate(zip(combo, self.parent_aa)):

            if mut != wt:
                mut_list.append(mut)
                # note the info dict positiosn is 1 indexed
                pos_list.append(self.lib_info["positions"][i+1])

        return mut_list, pos_list

    @property
    def df(self) -> pd.DataFrame:
            
        """
        Get the dataframe with mutation details
        """

        df = self.input_df.copy()
        df[[self._mut_col_name, self._pos_col_name]] = df[self.combo_col_name].apply(
            lambda x: pd.Series(self._append_mut_dets(x)),
            axis=1
        )

        return df.copy()

    @property
    def max_n_mut(self) -> int:

        """
        Get the maximum number of mutations
        """

        return self.df["n_mut"].max()


class ESM(LibData):

    """
    A class for generating ESM scores
    """

    def __init__(self,
        input_csv: str,
        scale_fit: str,
        esm_model_name: str = "esm2_t33_650M_UR50D",
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        mut_col_name: str = "mut",
        pos_col_name: str = "pos",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        seq_dir: str = "data/seq",
        esm_dir: str = "zs/esm",
        regen_esm=False
        ):

        super().__init__(input_csv, scale_fit, combo_col_name, mut_col_name, pos_col_name, seq_col_name, fit_col_name, seq_dir)

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

        self._esm_dir = checkNgen_folder(esm_dir)
        self._esm_output_dir = checkNgen_folder(os.path.join(self._esm_dir, "output"))
        self._esm_path = os.path.join(self._esm_output_dir, f"{self.lib_name}.csv")

        self._logit_dir = checkNgen_folder(os.path.join(self._esm_dir, "logits", self._esm_model_name))
        self._logit_path = os.path.join(self._logit_dir, f"{self.lib_name}.npy")

        if self._logits_path != "" and os.path.exists(self._logits_path) and not(regen_esm):
            print(f"{self._logits_path} exists and regen_esm = {regen_esm}. Loading...")
            self._logits = np.load(self._logits_path)
        else:
            print(f"Generating {self._logits_path}...")
            self._logits = self._get_logits()

        # generate esm score
        self._esm_df = self._get_esm_score()

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
                token_probs = torch.log_softmax(
                    self._model(batch_tokens_masked)["logits"], dim=-1
                ).cpu().numpy()

            logits[i] = token_probs[0, i+1]

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

    def _get_esm_score(self, df: pd.DataFrame, _sum: bool = False):

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

        if _sum:
            for i, mut in enumerate(df["mut"]):
                s = np.zeros(len(mut))
                for j, mt in enumerate(mut):
                    if mt == "WT":
                        score[i] = 0
                        continue

                    elif mt == "NA":
                        score[i] = np.nan
                        continue

                    else:
                        pos = (
                            int(df["pos"].iloc[i][j]) - 1
                        )  # Position of the mutation with python indexing
                        wt = wt_sequence[pos]
                        s[j] = self._get_mutant_prob(mt=mt, wt=wt, pos=pos)
                    
                    score[i] += s.sum()

        else:
            for i, mut in enumerate(df["mut"]):

                mt = mut[0]

                if mt == "WT":
                    score[i] = 0
                    continue

                elif mt == "NA":
                    score[i] = np.nan
                    continue

                else:
                    pos = int(df["pos"].iloc[i][0] - 1)
                    wt = wt_sequence[pos]
                    score[i] = self._get_mutant_prob(mt=mt, wt=wt, pos=pos)

        return score

    def _get_n_df(self, n: int = 1):

        """
        Get n data frame with n mutants
        """

        return self.df[self.df["n_mut"] == n].copy()

    def _get_esm_score(self):

        """
        Get any score for each variant in the data set
        """

        df_n_list = []

        # Get the n mutant scores
        for i in list(range(self.max_n_mut+1))[1:]:
            # Filter out n mutants
            df_n = self._get_n_df(i)
            if df_n.empty:  # Check if the DataFrame is empty after filtering
                assert "Data set is empty"
                continue

            if i == 1:
                score_n = self._get_esm_score(df_n, _sum=False)
            else:
                score_n = self._get_esm_score(df_n, _sum=True)

            # Add column with number of mutations
            df_n.loc[:, "esm_score"] = score_n

            df_n_list.append(df_n)

        return pd.concat(df_n_list, axis=0)

