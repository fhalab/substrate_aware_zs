"""
Test the chai processing.
"""

import sys
import os


from datetime import datetime

from REVIVAL.zs.af3 import run_af3_struct, parse_all_af3_scores
from REVIVAL.util import checkNgen_folder

if __name__ == "__main__":

    # os.environ["CUDA_VISIBLE_DEVICES"] = "1"

    # log outputs
    f = open(
        os.path.join(
            checkNgen_folder("logs/zs/af3"),
            f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out",
        ),
        "w",
    )
    sys.stdout = f

    # run_af3_struct(
    #     # "data/meta/not_scaled/Rma-*.csv",
    #     # "data/meta/not_scaled/PfTrpB-*.csv",
    #     [ 
    #     # #     # "data/meta/not_scaled/PfTrpB-4bromo.csv", 
    #     # #     # "data/meta/not_scaled/PfTrpB-4cyano.csv",
    #     # #     # "data/meta/not_scaled/PfTrpB-5bromo.csv",
    #     # #     # "data/meta/not_scaled/PfTrpB-5chloro.csv",
    #     # #     # "data/meta/not_scaled/PfTrpB-5cyano.csv",
    #     # #     # "data/meta/not_scaled/PfTrpB-5iodo.csv",
    #     # #     # "data/meta/not_scaled/PfTrpB-6chloro.csv",
    #     # #     # "data/meta/not_scaled/PfTrpB-7bromo.csv",
    #     # #     # "data/meta/not_scaled/PfTrpB-7iodo.csv",
    #     # #     # "data/meta/not_scaled/PfTrpB-7methyl.csv",
    #     # #     # "data/meta/not_scaled/PfTrpB-56chloro.csv",
    #         "data/meta/not_scaled/Rma-CB.csv",
    #         "data/meta/not_scaled/Rma-CSi.csv",
    #         "data/meta/not_scaled/ParLQ.csv",
    #     ],
    #     # gen_opt="joint",
    #     gen_opt="seperate",
    #     cofactor_dets="cofactor",
    #     gpu_id="0"
    # )

    # run_af3_struct(
    #     "data/meta/not_scaled/*.csv",
    #     gen_opt="joint",
    #     cofactor_dets="cofactor",
    #     gpu_id="0"
    # )

    parse_all_af3_scores(af3_struct_dir = "zs/af3/struct_joint")
    parse_all_af3_scores(af3_struct_dir = "zs/af3/struct_seperate")

    # def run_af3_struct(
    #     pattern: str | list = "data/meta/not_scaled/*.csv", 
    #     gen_opt: str = "joint",
    #     kwargs: dict = {}
    # ):
    # parse_all_af3_scores(af3_struct_dir: str = "zs/af3/struct_joint")

    f.close()