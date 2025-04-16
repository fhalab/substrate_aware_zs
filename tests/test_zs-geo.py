"""
Test extracting distances etc from docked chai or af3 structures.
"""

import sys
import os

from datetime import datetime

from substrate_aware.zs.geometries import run_bonddist
from substrate_aware.util import checkNgen_folder

if __name__ == "__main__":

    # log outputs
    f = open(
        os.path.join(
            checkNgen_folder("logs/zs/geo"),
            f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out",
        ),
        "w",
    )
    sys.stdout = f


    # run_bonddist(
    #     struct_dir="zs/af3/struct_joint",
    #     pattern="data/meta/PfTrpB*.csv",
    # )
    
    # run_bonddist(
    #     struct_dir="zs/af3/struct_seperate",
    #     pattern="data/meta/ParLQ-*.csv",
    #     )

    # run_bonddist(
    #     struct_dir="zs/af3/struct_seperate",
    #     pattern="data/meta/Rma-*.csv",
    #     )

    # run_bonddist(
    #     struct_dir="zs/chai/struct_joint",
    #     pattern="data/meta/PfTrpB*.csv",
    # )

    # run_bonddist(
    #     struct_dir="zs/chai/struct_seperate",
    #     pattern="data/meta/ParLQ-*.csv",
    #     )
    
    # run_bonddist(
    #     struct_dir="zs/chai/struct_seperate",
    #     pattern="data/meta/Rma-*.csv",
    #     )


    f.close()

    """
    run_bonddist(
        struct_dir: str,
        pattern: str | list = "data/meta/*.csv",
        kwargs: dict = {},
    )
    """