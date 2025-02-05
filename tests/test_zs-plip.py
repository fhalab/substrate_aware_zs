"""
Test plip
"""

import sys
import os
from glob import glob

from datetime import datetime

from REVIVAL.zs.plip import run_lib_plip, run_all_plip_zs
from REVIVAL.util import checkNgen_folder

if __name__ == "__main__":

    # log outputs
    f = open(
        os.path.join(
            checkNgen_folder("logs/zs/plip"),
            f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out",
        ),
        "w",
    )
    sys.stdout = f

    # run_lib_plip(in_dir="data/structure", out_dir="zs/plip")

    # for lib_dir in sorted(glob("zs/af3/struct_joint/*_scope")):
    #     run_lib_plip(in_dir=lib_dir, out_dir="zs/plip")

    for lib_dir in sorted(glob("zs/af3/struct_seperate/ParLQ-*")):
        run_lib_plip(in_dir=lib_dir, out_dir="zs/plip")

    # for lib_dir in sorted(glob("zs/chai/struct_joint/PfTrpB_scope")):
    #     run_lib_plip(in_dir=lib_dir, out_dir="zs/plip")

    # for lib_dir in sorted(glob("zs/chai/struct_seperate/PfTrpB_scope")):
    #     run_lib_plip(in_dir=lib_dir, out_dir="zs/plip")

    for plip_dir in[
        # "zs/plip/af3/struct_joint",
        # "zs/plip/chai/struct_joint",
        "zs/plip/af3/struct_seperate",
        # "zs/plip/chai/struct_seperate",
    ]:
    
        run_all_plip_zs(
            "data/meta/ParLQ-*.csv",
            plip_dir
            )

    
    # for plip_dir in[
    #     "zs/plip/af3/struct_seperate",
    #     # "zs/plip/chai/struct_seperate",
    # ]:
    
    #     run_all_plip_zs(
    #         [
    #             # "data/meta/ParLQ.csv",
    #             "data/meta/Rma-CB_scope.csv",
    #             "data/meta/Rma-CSi_scope.csv",
    #         ],
    #         plip_dir
    #         )

    f.close()

    """
    run_lib_plip(in_dir: str, out_dir: str="zs/plip"):

    
    Get plip report for each of the variant in a given directory 

    if  in_dir = 'data/structure'
        out_dir = 'zs/plip'
        will look for structure directly under the folder, i.e.
            data/structure/PfTrpB.pdb
            data/structure/Rma.pdb
        to generate plip results under holo subdirectory in the out_dir, i.e.
            zs/plip/holo/PfTrpB/
            zs/plip/holo/Rma/

    if  in_dir = zs/af3/struct_joint/ParLQ
        out_dir = zs/plip
        will look for structures under the subfolders for each variant, i.e.
            zs/af3/struct_joint/ParLQ/w56e_y57k_l59f_q60d_f89w/w56e_y57k_l59f_q60d_f89w_model.cif
            zs/af3/struct_joint/ParLQ/w56e_y57k_l59f_q60d_f89w/seed-1_sample-0/model.cif
            zs/af3/struct_joint/ParLQ/w56e_y57k_l59f_q60d_f89w/seed-1_sample-1/model.cif
        to first convert cif to pdb and then
        to generate plip results under the out_dir that 
        perserve the structure details as well as consolidate and rename the variants and reps, i.e.
            zs/plip/af3/struct_joint/ParLQ/W56E:Y57K:L59F:Q60D:F89W_agg/
            zs/plip/af3/struct_joint/ParLQ/W56E:Y57K:L59F:Q60D:F89W_0/
            zs/plip/af3/struct_joint/ParLQ/W56E:Y57K:L59F:Q60D:F89W_1/

    if  in_dir = zs/chai/struct_joint/ParLQ
        out_dir = zs/plip
        will look for structures under the subfolders for each variant, i.e.
            zs/chai/struct_joint/ParLQ/W56A:Y57C:L59S:Q60E:F89G/W56A:Y57C:L59S:Q60E:F89G_0.cif
            zs/chai/struct_joint/ParLQ/W56A:Y57C:L59S:Q60E:F89G/W56A:Y57C:L59S:Q60E:F89G_1.cif
        to first convert cif to pdb and then
        to generate plip results under the out_dir that
        perserve the structure details, i.e.
            zs/plip/chai/struct_joint/ParLQ/W56A:Y57C:L59S:Q60E:F89G_0/
            zs/plip/chai/struct_joint/ParLQ/W56A:Y57C:L59S:Q60E:F89G_1/

    run_all_plip_zs(pattern: str | list, plip_dir: str, kwargs: dict = {})
    """