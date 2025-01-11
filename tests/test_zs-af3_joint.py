"""
Test the chai processing.
"""

import sys
import os


from datetime import datetime

from REVIVAL.zs.af3 import run_af3_prep
from REVIVAL.util import checkNgen_folder

if __name__ == "__main__":

    # os.environ["CUDA_VISIBLE_DEVICES"] = "0"

    # log outputs
    f = open(
        os.path.join(
            checkNgen_folder("logs/zs/af3"),
            f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out",
        ),
        "w",
    )
    sys.stdout = f

    # run_gen_chai_structure("data/meta/not_scaled/*.csv")
    run_af3_prep(
        # "data/meta/not_scaled/PfTrpB-*.csv",
        [ 
        #     # "data/meta/not_scaled/PfTrpB-4bromo.csv", 
        #     # "data/meta/not_scaled/PfTrpB-4cyano.csv",
        #     # "data/meta/not_scaled/PfTrpB-5bromo.csv",
        #     # "data/meta/not_scaled/PfTrpB-5chloro.csv",
        #     # "data/meta/not_scaled/PfTrpB-5cyano.csv",
        #     # "data/meta/not_scaled/PfTrpB-5iodo.csv",
        #     # "data/meta/not_scaled/PfTrpB-6chloro.csv",
        #     # "data/meta/not_scaled/PfTrpB-7bromo.csv",
        #     # "data/meta/not_scaled/PfTrpB-7iodo.csv",
        #     # "data/meta/not_scaled/PfTrpB-7methyl.csv",
        #     # "data/meta/not_scaled/PfTrpB-56chloro.csv",
            # "data/meta/not_scaled/Rma-CB.csv",
            # "data/meta/not_scaled/Rma-CSi.csv",
            "data/meta/not_scaled/ParLQ.csv",
        ],
        gen_opt="joint",
        # gen_opt="seperate",
        cofactor_dets="cofactor",
        gpu_id="1"
    )

    # def run_af3_prep(
    #     pattern: str | list = "data/meta/not_scaled/*.csv", 
    #     gen_opt: str = "joint",
    #     kwargs: dict = {}
    # ):

    f.close()