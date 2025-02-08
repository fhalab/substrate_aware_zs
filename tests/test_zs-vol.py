"""
Test plip
"""

import sys

import os
from glob import glob

from datetime import datetime

from REVIVAL.zs.vol import run_all_vol
from REVIVAL.util import checkNgen_folder

if __name__ == "__main__":

    # log outputs
    f = open(
        os.path.join(
            checkNgen_folder("logs/zs/vol"),
            f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out",
        ),
        "w",
    )
    sys.stdout = f

   
    run_all_vol(pattern="data/meta/*_scope.csv", samesub=False)
    # run_all_vol(pattern="data/meta/PfTrpB-*.csv", samesub=True)
    # run_all_vol(pattern="data/meta/ParLQ-*.csv", samesub=True)
    # run_all_vol(pattern=
    #     ["data/meta/Rma-CB.csv", "data/meta/Rma-CSi.csv"], samesub=True)

    f.close()