"""
Test ev mutation ZS
"""

import sys
import os

from glob import glob

from datetime import datetime

from substrate_aware.zs.ev import run_all_ev
from substrate_aware.util import checkNgen_folder

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

    # run_all_ev("data/meta/*.csv")

    f.close()