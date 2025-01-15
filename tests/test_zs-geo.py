"""
Test the vina docking and extracting score.
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


    run_bonddist(
        struct_dir="zs/af3/struct_joint",
        pattern="data/meta/not_scaled/PfTrpB*.csv",
    )
    
    # run_bonddist(
    #     struct_dir="zs/af3/struct_seperate",
    #     pattern=[
    #         "data/meta/not_scaled/ParLQ.csv",
    #         "data/meta/not_scaled/Rma-CB.csv",
    #         "data/meta/not_scaled/Rma-CSi.csv",
    #         ]
    #     )

    run_bonddist(
        struct_dir="zs/chai/structure_joint",
        pattern="data/meta/not_scaled/PfTrpB*.csv",
    )

    # run_bonddist(
    #     struct_dir="zs/chai/structure_seperate",
    #     pattern=[
    #         "data/meta/not_scaled/ParLQ.csv",
    #         "data/meta/not_scaled/Rma-CB.csv",
    #         "data/meta/not_scaled/Rma-CSi.csv",
    #         ]
    #     )


    f.close()

    """
    run_bonddist(
        struct_dir: str,
        pattern: str | list = "data/meta/not_scaled/*.csv",
        kwargs: dict = {},
    )
    """