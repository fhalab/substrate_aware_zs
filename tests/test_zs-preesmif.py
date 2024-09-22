"""Test the triad pre and post processing."""

import sys
import os

from glob import glob

from datetime import datetime

from REVIVAL.zs.esmif import run_esmif_input_file
from REVIVAL.util import checkNgen_folder

if __name__ == "__main__":

    # log outputs
    f = open(
        os.path.join(
            checkNgen_folder("logs/zs/esmif/pre"),
            f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out",
        ),
        "w",
    )
    sys.stdout = f

    run_esmif_input_file("data/meta/scale2parent/*.csv")

    """run_esmif_input_file(
    pattern: str | list = "data/meta/scale2parent/*.csv"
    )"""

    f.close()