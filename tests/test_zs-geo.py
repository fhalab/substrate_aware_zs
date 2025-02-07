"""
Test extracting distances etc from docked chai or af3 structures.
"""

import sys
import os

from datetime import datetime

from REVIVAL.zs.geometries import run_bonddist
from REVIVAL.util import checkNgen_folder

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
    #     pattern="data/meta/*.csv",
    # )
    
    run_bonddist(
        struct_dir="zs/af3/struct_seperate",
        pattern="data/meta/ParLQ-*.csv",
        # pattern=[
        #     "data/meta/ParLQ.csv",
        #     "data/meta/Rma-CB.csv",
        #     "data/meta/Rma-CSi.csv",
        #     ]
        )

    # run_bonddist(
    #     struct_dir="zs/chai/struct_joint",
    #     pattern="data/meta/*.csv",
    # )

    run_bonddist(
        struct_dir="zs/chai/struct_seperate",
        pattern="data/meta/ParLQ-*.csv",
        # pattern=[
        #     "data/meta/ParLQ.csv",
        #     "data/meta/Rma-CB.csv",
        #     "data/meta/Rma-CSi.csv",
        #     ]
        )


    f.close()

    """
    run_bonddist(
        struct_dir: str,
        pattern: str | list = "data/meta/*.csv",
        kwargs: dict = {},
    )
    """