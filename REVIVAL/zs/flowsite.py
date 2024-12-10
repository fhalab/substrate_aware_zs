# File to automatically run Flowsite as a ZS
# Author: Lukas Radtke {lradtke@caltech.edu}
# Date: 06/12/24

import numpy as np
import pandas as pd
import os 
import subprocess
import re
import ipdb
import glob
from copy import deepcopy
from tqdm import tqdm
from scipy import stats

from REVIVAL.preprocess import ZSData
from REVIVAL.util import checkNgen_folder


class FlowsiteData(ZSData):
    def __init__(
            self,
            input_csv = str,
            scale_fit: str = 'parent',
            structure_dir: str = 'data/structure',
            flowsite_dir: str = '../FlowSite',
            results_dir: str = './outputs',
            combo_col_name: str = 'AAs',
            var_col_name: str = 'var',
            seq_col_name: str = 'seq',
            fit_col_name: str = 'fitness',
            substrate_col_name: str = 'substrate',
            cofactor_col_name: str = None,

    ):
        
        super.__init__(
            input_csv = input_csv,
            scale_fit = scale_fit,
            combo_col_name = combo_col_name,
            var_col_name = var_col_name,
            seq_col_name = seq_col_name,
            fit_col_name = fit_col_name,
            structure_dir = structure_dir,
            substrate_col_name = substrate_col_name,
            cofactor_col_name = cofactor_col_name,
        )

        self.structure_dir = structure_dir
        self.flowsite_dir = flowsite_dir   
        self.results_dir = results_dir     


def inference(self) -> None:
    """
    Prepares the inputs for the model and runs FlowSite
    """
    csv = pd.read_csv(self.input_csv, keep_default_na=False)
    csv = csv[~csv.apply(lambda row: row.astype(str).str.contains(r"\*").any(), axis=1)].reset_index(drop=True)  # filters for stop codons, marked as '*'

    if self.cofactor_col_name is not None:
        #generate concatenated smiles
        pass
    else:
        # only use single substrate smiles
        pass

def run_flowsite(pattern: str | list = None, kwargs: dict = {}):

    if isinstance(pattern, str):
        lib_list = glob(pattern)
    else:
        lib_list = deepcopy(pattern)

    for lib in lib_list:
        print(f'Running Flowsite for {lib}...')
        FlowsiteData(input_csv=lib, **kwargs)


# working: CUDA_VISIBLE_DEVICES="" python -m inference --num_inference 10 --out_dir data/output --design_residues "A165, A183, A301" --smiles "O=Cc1c(O)c(C)ncc1COP(=O)(O)O, "  --protein ../enzymeoracle/data/structure/PfTrpB-4bromo.pdb --pocket_def_residues "A165, A183, A301" --batch_size 16 --checkpoint pocket_gen/lf5t55w4/checkpoints/best.ckpt --run_test --run_name inference1 --layer_norm --fake_constant_dur 100000 --fake_decay_dur 10 --fake_ratio_start 0.2 --fake_ratio_end 0.2 --residue_loss_weight 0.2 --use_tfn --use_inv --time_condition_tfn --correct_time_condition --time_condition_inv --time_condition_repeat --flow_matching --flow_matching_sigma 0.5 --prior_scale 1 --num_angle_pred 5 --ns 48 --nv 12 --batch_norm --self_condition_x --self_condition_inv --no_tfn_self_condition_inv --self_fancy_init --pocket_residue_cutoff 12 
