"""
Test the esmif pre and post processing.
HAVE TO use esmif.yml to run this.
"""

import sys
import os

from datetime import datetime

from substrate_aware.zs.ligandmpnn import run_ligandmpnn, LigandMPNN_MODEL_DICT
from substrate_aware.util import checkNgen_folder

if __name__ == "__main__":

    os.environ["CUDA_VISIBLE_DEVICES"] = "1"

    # log outputs
    f = open(
        os.path.join(
            checkNgen_folder("logs/zs/ligandmpnn"),
            f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out",
        ),
        "w",
    )
    sys.stdout = f

    # for noise in LigandMPNN_MODEL_DICT.keys():
    #     run_ligandmpnn(
    #         pattern="data/meta/*.csv",
    #         noise_level=noise,
    #     )


    f.close()
