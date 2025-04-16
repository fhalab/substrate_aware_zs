"""Test the esm module."""

import sys
import os

from datetime import datetime

from substrate_aware.zs.esm import run_all_esm
from substrate_aware.util import checkNgen_folder

if __name__ == "__main__":
    
    os.environ["CUDA_VISIBLE_DEVICES"] = "1"

    log_folder = checkNgen_folder("logs/zs/esm")

    # log outputs
    f = open(os.path.join(log_folder, f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out"), 'w')
    sys.stdout = f

    # run_all_esm(pattern="data/meta/*")

    f.close()


"""
run_all_esm(
    pattern: str = "data/meta/*",
    scale_fit: str = "",
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
)
"""