"""
Test the chai processing.
"""

import sys
import os


from datetime import datetime

from REVIVAL.zs.chai import run_gen_chai_structure, parse_all_chai_scores
from REVIVAL.util import checkNgen_folder

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

    run_gen_chai_structure(
        "data/meta/ParLQ-*.csv",
        gen_opt="seperate",
        cofactor_dets="cofactor",
        samesub=True
        )

    # # run_gen_chai_structure("data/meta/*.csv")
    # run_gen_chai_structure([ 
    #     # "data/meta/PfTrpB-4bromo.csv",
    #     # "data/meta/PfTrpB-4bromo.csv", 
    #     # "data/meta/PfTrpB-4cyano.csv",
    #     # "data/meta/PfTrpB-5bromo.csv",
    #     # "data/meta/PfTrpB-5chloro.csv",
    #     # "data/meta/PfTrpB-5cyano.csv",
    #     # "data/meta/PfTrpB-5iodo.csv",
    #     # "data/meta/PfTrpB-6chloro.csv",
    #     # "data/meta/PfTrpB-7bromo.csv",
    #     # "data/meta/PfTrpB-7iodo.csv",
    #     # "data/meta/PfTrpB-7methyl.csv",
    #     # "data/meta/PfTrpB-56chloro.csv",
        

    #     "data/meta/Rma-CB.csv",
    #     "data/meta/Rma-CSi.csv",
    #     "data/meta/ParLQ.csv",
    # ], 
    # # gen_opt="joint-cofactor-no-substrate",
    # # cofactor_dets="transiminated-cofactor",
    # gen_opt="joint",
    # cofactor_dets="cofactor"
    # )

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
    #     "/disk2/fli/REVIVAL2/data/meta/Rma-CB_scope.csv",
    #     "/disk2/fli/REVIVAL2/data/meta/Rma-CSi_scope.csv",
    #     "/disk2/fli/REVIVAL2/data/meta/PfTrpB_scope.csv"
    # ],
    # gen_opt="joint",
    # cofactor_dets="cofactor",
    # samesub=False,
    # )

    # run_gen_chai_structure([
    #     "/disk2/fli/REVIVAL2/data/meta/Rma-CB_scope.csv",
    #     "/disk2/fli/REVIVAL2/data/meta/Rma-CSi_scope.csv",
    #     "/disk2/fli/REVIVAL2/data/meta/PfTrpB_scope.csv"
    # ],
    # gen_opt="seperate",
    # cofactor_dets="cofactor",
    # samesub=False,
    # )

    # run_gen_chai_structure([
    #     "/disk2/fli/REVIVAL2/data/substrate_scope/Rma-CSi.csv",
    #     # "/disk2/fli/REVIVAL2/data/substrate_scope/TrpB-platform.csv"
    # ],
    # gen_opt="joint",
    # cofactor_dets="cofactor",
    # samesub=False,
    # kwargs={"chai_dir": "chai_substratescope", "ifrerun": True}
    # )
    # run_gen_chai_structure([
    #     "/disk2/fli/REVIVAL2/data/substrate_scope/Rma-CSi.csv"
    #     # "/disk2/fli/REVIVAL2/data/substrate_scope/TrpB-platform.csv"
    # ],
    # gen_opt="seperate",
    # cofactor_dets="cofactor",
    # samesub=False,
    # kwargs={"chai_dir": "chai_substratescope", "ifrerun": True}
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