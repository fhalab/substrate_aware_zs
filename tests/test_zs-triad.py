"""Test the triad pre and post processing."""

import sys
import os

from glob import glob

from datetime import datetime

from substrate_aware.zs.triad import run_traid
from substrate_aware.util import checkNgen_folder

if __name__ == "__main__":

    # log outputs
    f = open(
        os.path.join(
            checkNgen_folder("logs/zs/triad"),
            f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out",
        ),
        "w",
    )
    sys.stdout = f

    # run_traid(
    #     "data/meta/*.csv",
    #     in_structure_dir="data/structure/docked",
    #     kwargs={"withsub": True, "cleanup": False}
    # )

    # run_traid("data/meta/*.csv", kwargs={"withsub": False, "cleanup": True})
    
    # run_traid(
    #     "data/meta/*.csv",
    #     in_structure_dir="data/structure/docked",
    #     kwargs={"withsub": True, "cleanup": True}
    # )


    f.close()