"""
Test the chai processing.
"""

import sys
import os


from datetime import datetime

from REVIVAL.zs.chai import run_gen_chai_structure, parse_all_chai_scores
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

    # # run_gen_chai_structure("data/meta/not_scaled/*.csv")
    # run_gen_chai_structure([ 
    #     # "data/meta/not_scaled/PfTrpB-4bromo.csv",
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
        

    #     "data/meta/not_scaled/Rma-CB.csv",
    #     "data/meta/not_scaled/Rma-CSi.csv",
    #     "data/meta/not_scaled/ParLQ.csv",
    # ], 
    # # gen_opt="joint-cofactor-no-substrate",
    # # cofactor_dets="transiminated-cofactor",
    # gen_opt="joint",
    # cofactor_dets="cofactor"
    # )

    parse_all_chai_scores(chai_struct_dir= "zs/chai/structure_joint")
    parse_all_chai_scores(chai_struct_dir= "zs/chai/structure_seperate")

    f.close()