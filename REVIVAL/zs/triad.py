"""
A script for preprocessing the data for the triad and extract energies
"""

import re
import os

from copy import deepcopy
from glob import glob

import pandas as pd
import numpy as np
from itertools import product

from REVIVAL.global_param import LIB_INFO_DICT
from REVIVAL.preprocess import ZSData
from REVIVAL.util import checkNgen_folder, get_file_name


class TriadData(ZSData):
    def __init__(self,
        input_csv: str,
        scale_fit: str,
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        mut_col_name: str = "mut",
        pos_col_name: str = "pos",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        seq_dir: str = "data/seq",
        zs_dir: str = "zs",
        triad_dir: str = "triad",
        chain_id: str = "A",
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

        self._triad_dir = checkNgen_folder(triad_dir)
        self._chain_id = chain_id


        
