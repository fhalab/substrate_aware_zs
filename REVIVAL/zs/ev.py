"""A script for running the EV mutation ZS data."""

from __future__ import annotations

import os
from itertools import chain

from glob import glob
from copy import deepcopy

from evcouplings.couplings import CouplingsModel

from REVIVAL.preprocess import ZSData
from REVIVAL.util import checkNgen_folder


class EVData(ZSData):
    """
    A class for generating EV mutation scores
    """

    def __init__(
        self,
        input_csv: str,
        scale_fit: str = "parent",
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

        # create the subfolder
        self._ev_dir = checkNgen_folder(os.path.join(self._zs_dir, ev_dir))

        print(f"EV data will be saved at {self._ev_dir}")

        # get the ev score
        self._ev_df = self._get_evscore()

        # save the ev score
        self._ev_df.to_csv(self.ev_csv, index=False)


    def _check_idx(self):
        """Check if the position is in the index map"""
        all_pos = set(chain.from_iterable(self.df[self._pos_col_name]))
        return all([p in self._idx_map for p in all_pos])

    def _get_evscore(self):

        df = deepcopy(self.df)
        
        def parse(x):
            # (self.target_positions[j], self.wt_aas[j], mut_char)
            return (int(x[1:-1]), x[0], x[-1])

        df['ev_score'] = df[self._var_col_name].apply(
            lambda x: 0 if x == "WT" else self._model.delta_hamiltonian([parse(m) for m in x.split(':')])[0]
        )
        
        return df

    @property
    def ev_csv(self) -> str:
        """
        A property for the triad csv
        """
        return os.path.join(self._ev_dir, f"{self.lib_name}.csv")

    
def run_all_ev(
    pattern: str | list = "data/meta/scale2parent/*.csv",
    kwargs: dict = {}
    ):

    """
    Run the ev mutation scores for all the libraries
    
    Args:
    """
    if isinstance(pattern, str):
        lib_list = glob(pattern)
    else:
        lib_list = deepcopy(pattern)

    for lib in lib_list:
        print(f"Getting ev zs for {lib}...")
        EVData(input_csv=lib, **kwargs)

