"""
Test the chai processing.
"""

import sys
import os


from datetime import datetime

from substrate_aware.zs.chai import run_gen_chai_structure, parse_all_chai_scores
from substrate_aware.util import checkNgen_folder

if __name__ == "__main__":

    os.environ["CUDA_VISIBLE_DEVICES"] = "0"

    # log outputs
    f = open(
        os.path.join(
            checkNgen_folder("logs/zs/chai/structure"),
            f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out",
        ),
        "w",
    )
    sys.stdout = f

    # run_gen_chai_structure(
    #     "data/meta/ParLQ-*.csv",
    #     gen_opt="seperate",
    #     cofactor_dets="cofactor",
    #     samesub=True
    #     )

    # run_gen_chai_structure(
    #     "data/meta/Rma-*.csv",
    #     gen_opt="seperate",
    #     cofactor_dets="cofactor",
    #     samesub=True
    #     )

    # run_gen_chai_structure(
    #     "data/meta/PfTrpB-*.csv",
    #     gen_opt="joint",
    #     cofactor_dets="cofactor",
    #     samesub=True
    #     )

    # run_gen_chai_structure([
    #     "data/meta/PfTrpB-4bromo.csv",
    #     "data/meta/Rma-CB.csv",
    #     "data/meta/Rma-CSi.csv",
    #     "data/meta/ParLQ.csv",
    # ],
    # gen_opt="apo"
    # )

    # run_gen_chai_structure([
    #     "data/meta/Rma-CB.csv",
    #     "data/meta/Rma-CSi.csv",
    #     "data/meta/ParLQ.csv",
    # ], 
    # gen_opt="joint-carbene_precursor-heme",
    # cofactor_dets="inactivated-cofactor",
    # samesub=True,
    # )

    # run_gen_chai_structure([
    #     "data/meta/Rma-CB.csv",
    #     "data/meta/Rma-CSi.csv",
    #     "data/meta/ParLQ.csv",
    # ], 
    # gen_opt="seperate-carbene_precursor-heme",
    # cofactor_dets="inactivated-cofactor",
    # samesub=True,
    # )
    
    # run_gen_chai_structure([
    #     # "data/meta/PfTrpB-4bromo.csv",
    #     # "data/meta/Rma-CB.csv",
    #     # "data/meta/Rma-CSi.csv",
    #     "data/meta/ParLQ.csv",
    # ],
    # gen_opt="joint-cofactor-no-substrate",
    # cofactor_dets="inactivated-cofactor"
    # )
    
    # parse_all_chai_scores(chai_struct_dir= "zs/chai/struct_joint")
    # parse_all_chai_scores(chai_struct_dir= "zs/chai/struct_seperate")

    f.close()