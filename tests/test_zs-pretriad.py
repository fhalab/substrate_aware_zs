"""Test the triad pre and post processing."""

import sys
import os

from glob import glob

from datetime import datetime

from REVIVAL.zs.triad import run_traid_gen_mut_file
from REVIVAL.util import checkNgen_folder

if __name__ == "__main__":

    # log outputs
    f = open(
        os.path.join(
            checkNgen_folder("logs/zs/triad/pre"),
            f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out",
        ),
        "w",
    )
    sys.stdout = f

    # run_traid_gen_mut_file("data/meta/not_scaled/D*.csv", kwargs={"withsub": False})
    run_traid_gen_mut_file("data/meta/not_scaled/TmT*.csv", kwargs={"withsub": False})

    """run_traid_gen_mut_file(
    pattern: str | list = "data/meta/scale2parent/*.csv"
    )"""

    f.close()