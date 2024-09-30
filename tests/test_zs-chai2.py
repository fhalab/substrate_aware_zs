"""
Test the chai processing.
"""

import sys
import os


from datetime import datetime

from REVIVAL.zs.chai import run_gen_chai_structure
from REVIVAL.util import checkNgen_folder

if __name__ == "__main__":

    os.environ["CUDA_VISIBLE_DEVICES"] = "1"

    # log outputs
    f = open(
        os.path.join(
            checkNgen_folder("logs/zs/chai/structure"),
            f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out",
        ),
        "w",
    )
    sys.stdout = f

    # run_gen_chai_structure("data/meta/scale2parent/*.csv")
    run_gen_chai_structure([
        # "data/meta/scale2parent/PfTrpB-5cyano.csv",
        # "data/meta/scale2parent/PfTrpB-56chloro.csv",
        # "data/meta/scale2parent/PfTrpB-7methyl.csv",
        "data/meta/scale2parent/PfTrpB-7iodo.csv",
    ])

    f.close()