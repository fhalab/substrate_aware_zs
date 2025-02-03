"""
Test the esmif pre and post processing.
HAVE TO use esmif.yml to run this.
"""

import sys
import os

from datetime import datetime

from REVIVAL.zs.flowsite import run_flowsite
from REVIVAL.util import checkNgen_folder

if __name__ == "__main__":

    os.environ["CUDA_VISIBLE_DEVICES"] = "1"

    # log outputs
    f = open(
        os.path.join(
            checkNgen_folder("logs/zs/flowsite"),
            f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out",
        ),
        "w",
    )
    sys.stdout = f

    for opt in [1,2]:
        run_flowsite(
            pattern="data/meta/*.csv",
            flowsite_inference_opt = "pocket_def_residues",
            flowsite_model_opt=opt
        )
    
    f.close()
