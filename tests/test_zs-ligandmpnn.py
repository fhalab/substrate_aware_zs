"""
Test the esmif pre and post processing.
HAVE TO use esmif.yml to run this.
"""

import sys
import os

from datetime import datetime

from REVIVAL.zs.ligandmpnn import run_ligandmpnn
from REVIVAL.util import checkNgen_folder

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

    #basedir is REVIVAL2
    run_ligandmpnn(
        pattern="data/meta/not_scaled/*.csv", 
        kwargs={'batch_size': 1}
        )

    """
    run_esmif(pattern: str | list = "data/meta/scale2parent/*.csv", kwargs: dict = {})
    """

    f.close()
