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

    for struct_dir in [
        "zs/chai/struct_joint/ParLQ",
        "zs/chai/struct_joint/Rma-CB",
        "zs/chai/struct_joint/Rma-CSi",
    ]:
        dock_lib_parallel(
            struct_dir = struct_dir,
            dock_opt="substrate+carbene",
            score_only=True,
            rerun=False,
            cofactor_dets="cofactor",
            max_workers=48
        )
    
    for struct_dir in [
        "zs/chai/struct_seperate/ParLQ",
        "zs/chai/struct_seperate/Rma-CB",
        "zs/chai/struct_seperate/Rma-CSi",
    ]:
        dock_lib_parallel(
            struct_dir = struct_dir,
            dock_opt="substrate+carbene",
            score_only=True,
            rerun=False,
            cofactor_dets="cofactor",
            max_workers=48
        )
    

    for struct_dir in [
        "zs/chai/struct_joint/ParLQ",
        "zs/chai/struct_joint/Rma-CB",
        "zs/chai/struct_joint/Rma-CSi",
    ]:
        dock_lib_parallel(
            struct_dir = struct_dir,
            dock_opt="substrate+carbene",
            score_only=False,
            rerun=False,
            cofactor_dets="cofactor",
            max_workers=48
        )
    
    for struct_dir in [
        "zs/chai/struct_seperate/ParLQ",
        "zs/chai/struct_seperate/Rma-CB",
        "zs/chai/struct_seperate/Rma-CSi",
    ]:
        dock_lib_parallel(
            struct_dir = struct_dir,
            dock_opt="substrate+carbene",
            score_only=False,
            rerun=False,
            cofactor_dets="cofactor",
            max_workers=48
        )
    


    for struct_dir in [
        "zs/af3/struct_joint/ParLQ",
        "zs/af3/struct_joint/Rma-CB",
        "zs/af3/struct_joint/Rma-CSi",
    ]:
        dock_lib_parallel(
            struct_dir = struct_dir,
            dock_opt="substrate+carbene",
            score_only=True,
            rerun=False,
            cofactor_dets="cofactor",
            max_workers=48
        )
    
    for struct_dir in [
        "zs/af3/struct_seperate/ParLQ",
        "zs/af3/struct_seperate/Rma-CB",
        "zs/af3/struct_seperate/Rma-CSi",
    ]:
        dock_lib_parallel(
            struct_dir = struct_dir,
            dock_opt="substrate+carbene",
            score_only=True,
            rerun=False,
            cofactor_dets="cofactor",
            max_workers=48
        )
    

    for struct_dir in [
        "zs/af3/struct_joint/ParLQ",
        "zs/af3/struct_joint/Rma-CB",
        "zs/af3/struct_joint/Rma-CSi",
    ]:
        dock_lib_parallel(
            struct_dir = struct_dir,
            dock_opt="substrate+carbene",
            score_only=False,
            rerun=False,
            cofactor_dets="cofactor",
            max_workers=48
        )
    
    for struct_dir in [
        "zs/af3/struct_seperate/ParLQ",
        "zs/af3/struct_seperate/Rma-CB",
        "zs/af3/struct_seperate/Rma-CSi",
    ]:
        dock_lib_parallel(
            struct_dir = struct_dir,
            dock_opt="substrate+carbene",
            score_only=False,
            rerun=False,
            cofactor_dets="cofactor",
            max_workers=48
        )
    
    f.close()
