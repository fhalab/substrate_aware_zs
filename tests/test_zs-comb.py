"""Test combing zs."""

import sys
import os

from datetime import datetime

from REVIVAL.zs.comb import run_all_combzs
from REVIVAL.util import checkNgen_folder

if __name__ == "__main__":
    
    os.environ["CUDA_VISIBLE_DEVICES"] = "1"

    log_folder = checkNgen_folder("logs/zs/comb")

    # log outputs
    f = open(os.path.join(log_folder, f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out"), 'w')
    sys.stdout = f

    # run_all_combzs(pattern="data/meta/*_scope.csv")
    # run_all_combzs(pattern="data/meta/*.csv")
    run_all_combzs(pattern="data/meta/ParLQ-*.csv")
        # [
        #     "data/meta/ParLQ.csv",
        #     "data/meta/PfTrpB-4bromo.csv",
        #     "data/meta/PfTrpB-4cyano.csv",
        #     "data/meta/PfTrpB-5bromo.csv",
        #     "data/meta/PfTrpB-5chloro.csv",
        #     "data/meta/PfTrpB-5cyano.csv",
        #     "data/meta/PfTrpB-5iodo.csv",
        #     "data/meta/PfTrpB-6chloro.csv",
        #     "data/meta/PfTrpB-7bromo.csv",
        #     "data/meta/PfTrpB-7iodo.csv",
        #     "data/meta/PfTrpB-7methyl.csv",
        #     "data/meta/PfTrpB-56chloro.csv",
        #     "data/meta/Rma-CB.csv",
        #     "data/meta/Rma-CSi.csv",
        # ]


    f.close()


"""
def run_all_combzs(
    pattern: str = "data/meta/*",
    scale_fit: str = "",
    combo_col_name: str = "AAs",
    var_col_name: str = "var",
    mut_col_name: str = "mut",
    pos_col_name: str = "pos",
    seq_col_name: str = "seq",
    fit_col_name: str = "fitness",
    seq_dir: str = "data/seq",
    zs_dir: str = "zs",
    comb_dir: str = "comb",
    zs_subdir_list: list = ["esm/output"],
)
"""