"""
Test the vina docking and extracting score.
"""

import sys
import os

from glob import glob

from datetime import datetime

from REVIVAL.zs.vina import dock_lib_parallel
from REVIVAL.util import checkNgen_folder

if __name__ == "__main__":

    # log outputs
    f = open(
        os.path.join(
            checkNgen_folder("logs/zs/vina/dock"),
            f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out",
        ),
        "w",
    )
    sys.stdout = f

    for struct_dir in sorted(glob(("zs/chai/struct_joint/PfTrpB-5*"))):
        dock_lib_parallel(
            struct_dir = struct_dir,
            dock_opt="substrate",
            score_only=True,
            rerun=False,
            cofactor_dets="cofactor",
        )
    
    for struct_dir in sorted(glob(("zs/chai/struct_joint/PfTrpB-6*"))):
        dock_lib_parallel(
            struct_dir = struct_dir,
            dock_opt="substrate",
            score_only=True,
            rerun=False,
            cofactor_dets="cofactor",
        )

    for struct_dir in sorted(glob(("zs/chai/struct_joint/PfTrpB-7*"))):
        dock_lib_parallel(
            struct_dir = struct_dir,
            dock_opt="substrate",
            score_only=True,
            rerun=False,
            cofactor_dets="cofactor",
        )

    for struct_dir in sorted(glob(("zs/chai/struct_seperate/PfTrpB-*"))):
        dock_lib_parallel(
            struct_dir = struct_dir,
            dock_opt="substrate",
            score_only=True,
            rerun=False,
            cofactor_dets="cofactor",
        )

    f.close()
