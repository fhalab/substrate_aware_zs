"""A script for running the EV mutation ZS data."""

from __future__ import annotations

import os
from itertools import chain

from glob import glob
from copy import deepcopy

from tqdm import tqdm

from evcouplings.couplings import CouplingsModel

from substrate_aware.preprocess import ZSData
from substrate_aware.util import checkNgen_folder, get_file_name


EV_META = {
    "ParLQ": {
        "recommended": {
            "bitscore": 0.1,
            "sequences": 15086,
            "seqs_per_l": 123.7,
            "quality": 10,
        },
        "chosen": {
            "bitscore": 0.7,
            "sequences": 343,
            "seqs_per_l": 1.8,
            "quality": 8,
        },
        "other_2": {
            "bitscore": 0.3,
            "sequences": 875,
            "seqs_per_l": 5.4,
            "quality": 6,
        },
        "other_3": {
            "bitscore": 0.5,
            "sequences": 343,
            "seqs_per_l": 1.8,
            "quality": 8,
        },
    },
    "PfTrpB": {
        "recommended": {
            "bitscore": 0.1,
            "sequences": 74795,
            "seqs_per_l": 262.4,
            "quality": 10,
        },
        "chosen": {
            "bitscore": 0.3,
            "sequences": 5996,
            "seqs_per_l": 15.7,
            "quality": 10,
        },
        "other_2": {
            "bitscore": 0.5,
            "sequences": 5935,
            "seqs_per_l": 15.6,
            "quality": 10,
        },
        "other_3": {
            "bitscore": 0.7,
            "sequences": 4647,
            "seqs_per_l": 12.2,
            "quality": 10,
        },
    },
    "Rma": {
        "recommended": {
            "bitscore": 0.3,
            "sequences": 79025,
            "seqs_per_l": 987.8,
            "quality": 10,
        },
        "chosen": {
            "bitscore": 0.5,
            "sequences": 3042,
            "seqs_per_l": 33.1,
            "quality": 10,
        },
        "other_3": {
            "bitscore": 0.7,
            "sequences": 1940,
            "seqs_per_l": 21.1,
            "quality": 10,
        },
    },
}


class EVData(ZSData):
    """
    A class for generating EV mutation scores
    """

    def __init__(
        self,
        input_csv: str,
        ev_model_dir: str = "data/evmodel",
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        mut_col_name: str = "mut",
        pos_col_name: str = "pos",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        protein_name: str = "",
        seq_dir: str = "data/seq",
        zs_dir: str = "zs",
        ev_dir: str = "ev",
    ):

        super().__init__(
            input_csv=input_csv,
            combo_col_name=combo_col_name,
            var_col_name=var_col_name,
            mut_col_name=mut_col_name,
            pos_col_name=pos_col_name,
            seq_col_name=seq_col_name,
            fit_col_name=fit_col_name,
            protein_name=protein_name,
            seq_dir=seq_dir,
            zs_dir=zs_dir,
        )

        self._ev_model_path = os.path.join(ev_model_dir, self.protein_name + ".model")
        print(f"Loading model at {self._ev_model_path}...")

        self._model = CouplingsModel(self._ev_model_path)
        print("Model loaded")
        self._idx_map = self._model.index_map

        # check if the position is in the index map
        # assert self._check_idx(), "Not all positions are in the index map!!!"

        # create the subfolder
        self._ev_dir = checkNgen_folder(os.path.join(zs_dir, ev_dir))

        print(f"EV data will be saved at {self._ev_dir}")

        # get the ev score
        self._ev_df = self._get_evscore()

        # save the ev score
        self._ev_df.to_csv(self.ev_csv, index=False)

    def _check_idx(self):
        """Check if the position is in the index map"""
        all_pos = set(chain.from_iterable(self.df[self._pos_col_name]))
        print(all_pos)
        print(self._idx_map)
        return all([p in self._idx_map for p in all_pos])

    def _get_evscore(self):

        df = deepcopy(self.df)

        def parse(x):
            # (self.target_positions[j], self.wt_aas[j], mut_char)
            return (int(x[1:-1]), x[0], x[-1])

        df["ev_score"] = df[self._var_col_name].apply(
            lambda x: 0
            if x == "WT"
            else self._model.delta_hamiltonian([parse(m) for m in x.split(":")])[0]
        )

        return df

    @property
    def ev_csv(self) -> str:
        """
        A property for the ev csv
        """
        return os.path.join(self._ev_dir, f"{self.lib_name}.csv")


def run_all_ev(pattern: str | list = "data/meta/*.csv", kwargs: dict = {}):

    """
    Run the ev mutation scores for all the libraries

    Args:
    """

    if isinstance(pattern, str):
        lib_list = glob(pattern)
    else:
        lib_list = deepcopy(pattern)

    for lib in tqdm(lib_list):
        print(f"Getting ev zs for {lib}...")
        EVData(input_csv=lib, **kwargs)