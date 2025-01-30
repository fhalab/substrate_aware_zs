"""
Test plip
"""

import sys

import os
from glob import glob

from datetime import datetime

from REVIVAL.zs.hydro import run_all_hydro
from REVIVAL.util import checkNgen_folder

if __name__ == "__main__":

    # log outputs
    f = open(
        os.path.join(
            checkNgen_folder("logs/zs/hydro"),
            f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out",
        ),
        "w",
    )
    sys.stdout = f

   
    # run_all_hydro(pattern="data/meta/not_scaled/*", plip_dir="zs/plip/af3/struct_joint")
    # run_all_hydro(pattern=[
    #     "data/meta/not_scaled/ParLQ.csv",
    #     "data/meta/not_scaled/Rma-CB.csv",
    #     "data/meta/not_scaled/Rma-CSi.csv",
    #     ], plip_dir="zs/plip/af3/struct_seperate")

    # run_all_hydro(pattern="data/meta/not_scaled/*", plip_dir="zs/plip/chai/struct_joint")
    # run_all_hydro(pattern=[
    #     "data/meta/not_scaled/ParLQ.csv",
    #     "data/meta/not_scaled/Rma-CB.csv",
    #     "data/meta/not_scaled/Rma-CSi.csv",
    #     ], plip_dir="zs/plip/chai/struct_seperate")

    run_all_hydro(pattern="data/meta/not_scaled/*", plip_dir="zs/plip/frompdb/struct")


    f.close()

    """
    run_hydro(pattern: str | list, plip_dir: str, kwargs: dict = {})
    """