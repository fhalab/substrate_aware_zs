"""
Test the triad pre and post processing.
HAVE TO use esmif.yml to run this.
"""

import sys
import os

from glob import glob

from datetime import datetime

from REVIVAL.zs.esmif import run_esmif
from REVIVAL.util import checkNgen_folder

if __name__ == "__main__":

    os.environ["CUDA_VISIBLE_DEVICES"] = "0"

    # log outputs
    f = open(
        os.path.join(
            checkNgen_folder("logs/zs/esmif"),
            f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out",
        ),
        "w",
    )
    sys.stdout = f

    # run_esmif("data/meta/scale2parent/5DW0*", kwargs={"structure_dir": "data/structure"})
    run_esmif("data/meta/not_scaled/amiE-i*", kwargs={"structure_dir": "data/structure", "withsub": False})

    """
    run_esmif(pattern: str | list = "data/meta/scale2parent/*.csv", kwargs: dict = {})
    """

    f.close()