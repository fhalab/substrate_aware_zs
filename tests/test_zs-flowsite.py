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

    run_flowsite(
        pattern="data/meta/not_scaled/*.csv",
    )
        
    f.close()
