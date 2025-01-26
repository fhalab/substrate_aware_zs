"""Test the triad pre and post processing."""

import sys
import os

from glob import glob

from datetime import datetime

from REVIVAL.zs.triad import run_traid
from REVIVAL.util import checkNgen_folder

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

    # run_traid("data/meta/not_scaled/*.csv", kwargs={"withsub": False, "cleanup": False})
    run_traid(
        "data/meta/not_scaled/*.csv",
        in_structure_dir="data/structure/docked",
        kwargs={"withsub": True, "cleanup": False}
    )


    f.close()