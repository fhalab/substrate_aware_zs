"""
Test ev mutation ZS
"""

import sys
import os

from glob import glob

from datetime import datetime

from REVIVAL.zs.ev import run_all_ev
from REVIVAL.util import checkNgen_folder

if __name__ == "__main__":

    os.environ["CUDA_VISIBLE_DEVICES"] = "0"

    # log outputs
    f = open(
        os.path.join(
            checkNgen_folder("logs/zs/ev"),
            f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out",
        ),
        "w",
    )
    sys.stdout = f

    # run_all_ev("data/meta/scale2parent/*.csv")
    kwargs = {
        "scale_fit": "max",
        "var_col_name": "muts",
        "ev_model_dir": "data/SSMuLA_ev",
        "zs_dir": "SSMuLA-zs",
        "ev_dir": "ev",
    }
    run_all_ev("data/SSMuLA/*", kwargs=kwargs)

    f.close()