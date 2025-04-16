"""Test the triad pre and post processing."""

import sys
import os

from glob import glob

from datetime import datetime

from substrate_aware.zs.triad import run_parse_triad_results
from substrate_aware.util import checkNgen_folder

if __name__ == "__main__":

    # log outputs
    f = open(
        os.path.join(
            checkNgen_folder("logs/zs/triad/post"),
            f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out",
        ),
        "w",
    )
    sys.stdout = f

    run_parse_triad_results("data/meta/*.csv")

    """run_parse_triad_results(
    pattern: str | list = "data/meta/scale2parent/*.csv"
    )"""

    f.close()