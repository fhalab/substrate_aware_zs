"""A script for running the EV mutation ZS data."""

from __future__ import annotations

import os
from glob import glob


from __future__ import annotations

import os
from glob import glob
from copy import deepcopy

import numpy as np
import pandas as pd

import torch

from evcouplings.couplings import CouplingsModel
from evcouplings.align.alignment import Alignment

from REVIVAL.preprocess import ZSData
from REVIVAL.util import checkNgen_folder


class EVData(ZSData):
    """
    A class for generating EV mutation scores
    """

    def __init__(
        self,
        input_csv: str,
        scale_fit: str,
        ev_model_dir: str = "data/evmodel",
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        mut_col_name: str = "mut",
        pos_col_name: str = "pos",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        seq_dir: str = "data/seq",
        zs_dir: str = "zs",
        ev_dir: str = "ev",
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
            zs_dir
        )

        self._ev_model_path = os.path.join(ev_model_dir, self.protein_name + ".model")
        print(f"Loading model at {self._ev_model_path}...")

        self._model = CouplingsModel(self._ev_model_path)
        print("Model loaded")
        self._idx_map = self._model.index_map

        # check if the position is in the index map
        assert self._check_idx(), "Not all positions are in the index map!!!"


    def _check_idx(self):
        """Check if the position is in the index map"""

        all_pos = self.df[self.pos_col_name].unique()
        return all([p in self.idx_map for p in all_pos])

    def _get_evscore(self):
        
        def parse(x):
            # (self.target_positions[j], self.wt_aas[j], mut_char)
            return (int(x[1:-1]), x[0], x[-1])

        self.df['ev_score'] = self.df["var_col_name"].apply(lambda x: self._model.delta_hamiltonian([parse(m) for m in x.split(':')])[0])
        
        return self.df
