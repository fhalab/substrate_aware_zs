"""
Test the chai processing.
"""

import sys
import os


from datetime import datetime

from REVIVAL.zs.af3 import run_af3_struct, parse_all_af3_scores
from REVIVAL.util import checkNgen_folder

if __name__ == "__main__":

    os.environ["CUDA_VISIBLE_DEVICES"] = "1"

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
    #     "data/meta/ParLQ-*.csv",
    #     gen_opt="seperate",
    #     cofactor_dets="cofactor",
    #     samesub=True,
    #     gpu_id="1"
    # )
    # run_af3_struct(
    #     # "data/meta/Rma-*.csv",
    #     # "data/meta/PfTrpB-*.csv",
    #     [ 
    #     # #     # "data/meta/PfTrpB-4bromo.csv", 
    #     # #     # "data/meta/PfTrpB-4cyano.csv",
    #     # #     # "data/meta/PfTrpB-5bromo.csv",
    #     # #     # "data/meta/PfTrpB-5chloro.csv",
    #     # #     # "data/meta/PfTrpB-5cyano.csv",
    #     # #     # "data/meta/PfTrpB-5iodo.csv",
    #     # #     # "data/meta/PfTrpB-6chloro.csv",
    #     # #     # "data/meta/PfTrpB-7bromo.csv",
    #     # #     # "data/meta/PfTrpB-7iodo.csv",
    #     # #     # "data/meta/PfTrpB-7methyl.csv",
    #     # #     # "data/meta/PfTrpB-56chloro.csv",
    #         "data/meta/Rma-CB.csv",
    #         "data/meta/Rma-CSi.csv",
    #         "data/meta/ParLQ.csv",
    #     ],
    #     # gen_opt="joint",
    #     gen_opt="seperate",
    #     cofactor_dets="cofactor",
    #     gpu_id="0"
    # )

    # run_af3_struct(
    #     "data/meta/*.csv",
    #     gen_opt="joint",
    #     cofactor_dets="cofactor",
    #     gpu_id="0"
    # )

    # run_af3_struct([
    #     "data/meta/PfTrpB-4bromo.csv",
    #     "data/meta/Rma-CB.csv",
    #     "data/meta/Rma-CSi.csv",
    #     "data/meta/ParLQ.csv",
    # ],
    # gen_opt="apo",
    # gpu_id="0"
    # )

    # run_af3_struct([
    #     "data/meta/PfTrpB-4bromo.csv",
    #     "data/meta/Rma-CB.csv",
    #     "data/meta/Rma-CSi.csv",
    #     "data/meta/ParLQ.csv",
    # ],
    # gen_opt="joint-cofactor-no-substrate",
    # cofactor_dets="inactivated-cofactor",
    # gpu_id="0"
    # )

    # run_af3_struct([
    #     "data/meta/Rma-CB.csv",
    #     "data/meta/Rma-CSi.csv",
    #     "data/meta/ParLQ.csv",
    # ],
    # gen_opt="joint-carbene_precursor-heme",
    # cofactor_dets="inactivated-cofactor",
    # samesub=True,
    # gpu_id="0"
    # )

    # run_af3_struct([
    #     "data/meta/Rma-CB.csv",
    #     "data/meta/Rma-CSi.csv",
    #     "data/meta/ParLQ.csv",
    # ],
    # gen_opt="seperate-carbene_precursor-heme",
    # cofactor_dets="inactivated-cofactor",
    # samesub=True,
    # gpu_id="0"
    # )

    # run_af3_struct([
    #     "/disk2/fli/REVIVAL2/data/meta/Rma-CB_scope.csv",
    #     "/disk2/fli/REVIVAL2/data/meta/Rma-CSi_scope.csv",
    #     "/disk2/fli/REVIVAL2/data/meta/PfTrpB_scope.csv"
    # ],
    # gen_opt="joint",
    # cofactor_dets="cofactor",
    # gpu_id="0",
    # samesub=False,
    # )

    # run_af3_struct([
    #     "/disk2/fli/REVIVAL2/data/meta/Rma-CB_scope.csv",
    #     "/disk2/fli/REVIVAL2/data/meta/Rma-CSi_scope.csv",
    # ],
    # gen_opt="seperate",
    # cofactor_dets="cofactor",
    # gpu_id="0",
    # samesub=False,
    # )

    # parse_all_af3_scores(af3_struct_dir = "zs/af3/struct_joint")
    # parse_all_af3_scores(af3_struct_dir = "zs/af3/struct_seperate")

    # def run_af3_struct(
    #     pattern: str | list = "data/meta/*.csv", 
    #     gen_opt: str = "joint",
    #     kwargs: dict = {}
    # ):
    # parse_all_af3_scores(af3_struct_dir: str = "zs/af3/struct_joint")
    parse_all_af3_scores(af3_struct_pattern = "zs/af3/struct_seperate/ParLQ-*")

    f.close()