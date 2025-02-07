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

   
    # run_all_hydro(pattern="data/meta/*", plip_dir="zs/plip/af3/struct_joint", kwargs={"max_workers": 4})
    # run_all_hydro(pattern=[
    #     "data/meta/ParLQ.csv",
    #     "data/meta/Rma-CB.csv",
    #     "data/meta/Rma-CSi.csv",
    #     ], plip_dir="zs/plip/af3/struct_seperate", kwargs={"max_workers": 4})

    # run_all_hydro(pattern="data/meta/*", plip_dir="zs/plip/chai/struct_joint", kwargs={"max_workers": 4})
    # run_all_hydro(pattern=[
    #     "data/meta/ParLQ.csv",
    #     "data/meta/Rma-CB.csv",
    #     "data/meta/Rma-CSi.csv",
    #     ], plip_dir="zs/plip/chai/struct_seperate", kwargs={"max_workers": 4})

    # run_all_hydro(pattern="data/meta/ParLQ.csv", plip_dir="zs/plip/frompdb/struct", kwargs={"max_workers": 4})
    
    # run_all_hydro(pattern="data/meta/*scope.csv", plip_dir="zs/plip/frompdb/struct", kwargs={"max_workers": 4})

    # run_all_hydro(pattern="data/meta/PfTrpB_scope.csv", plip_dir="zs/plip/af3/struct_joint", kwargs={"max_workers": 4})
    
    # run_all_hydro(pattern="data/meta/ParLQ-*.csv", plip_dir="zs/plip/frompdb/struct", kwargs={"max_workers": 4})
    run_all_hydro(pattern="data/meta/ParLQ-*.csv", plip_dir="zs/plip/af3/struct_seperate", kwargs={"max_workers": 4})



    f.close()

    """
    run_hydro(pattern: str | list, plip_dir: str, kwargs: dict = {})
    """