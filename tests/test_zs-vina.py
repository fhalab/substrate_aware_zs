"""
Test the vina docking and extracting score.
"""

import sys
import os

from glob import glob

from datetime import datetime

from substrate_aware.zs.vina import dock_lib_parallel, VinaLibDock
from substrate_aware.util import checkNgen_folder

if __name__ == "__main__":

    # os.environ["CUDA_VISIBLE_DEVICES"] = "1"

    # log outputs
    f = open(
        os.path.join(
            checkNgen_folder("logs/zs/vina/dock"),
            f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out",
        ),
        "w",
    )
    sys.stdout = f

    

    # for struct_dir in sorted(glob(("zs/af3/struct_joint/*"))):
    #     dock_lib_parallel(
    #         struct_dir = struct_dir,
    #         dock_opt="substrate",
    #         score_only=True,
    #         rerun=False,
    #         cofactor_dets="cofactor",
    #         max_workers=32
    #     )
    
    # for struct_dir in sorted(glob(("zs/af3/struct_seperate/*"))):
    #     dock_lib_parallel(
    #         struct_dir = struct_dir,
    #         dock_opt="substrate",
    #         score_only=True,
    #         rerun=False,
    #         cofactor_dets="cofactor",
    #         max_workers=32
    #     )

    # for struct_dir in sorted(glob(("zs/af3/struct_joint/*"))):
    #     dock_lib_parallel(
    #         struct_dir = struct_dir,
    #         dock_opt="substrate",
    #         score_only=False,
    #         rerun=False,
    #         cofactor_dets="cofactor",
    #         max_workers=32
    #     )
    
    # for struct_dir in sorted(glob(("zs/af3/struct_seperate/*"))):
    #     dock_lib_parallel(
    #         struct_dir = struct_dir,
    #         dock_opt="substrate",
    #         score_only=False,
    #         rerun=False,
    #         cofactor_dets="cofactor",
    #         max_workers=32
    #     )


    # for struct_dir in sorted(glob(("zs/af3/struct_seperate/ParLQ-*"))):
    #     dock_lib_parallel(
    #         struct_dir = struct_dir,
    #         dock_opt="substrate",
    #         score_only=True,
    #         rerun=False,
    #         cofactor_dets="cofactor",
    #         max_workers=16
    #     )

    
    # for struct_dir in sorted(glob(("zs/chai/struct_seperate/ParLQ-*"))):
    #     dock_lib_parallel(
    #         struct_dir = struct_dir,
    #         dock_opt="substrate",
    #         score_only=True,
    #         rerun=False,
    #         cofactor_dets="cofactor",
    #         max_workers=16
    #     )

    # for struct_dir in sorted(glob(("zs/chai/struct_joint/PfTrpB*"))):
    #     dock_lib_parallel(
    #         struct_dir = struct_dir,
    #         dock_opt="substrate",
    #         score_only=True,
    #         rerun=False,
    #         cofactor_dets="cofactor",
    #         max_workers=32
    #     )
    
    # for struct_dir in sorted(glob(("zs/chai/struct_seperate/ParLQ-*"))):
    #     dock_lib_parallel(
    #         struct_dir = struct_dir,
    #         dock_opt="substrate",
    #         score_only=True,
    #         rerun=False,
    #         cofactor_dets="cofactor",
    #         max_workers=32
    #     )

    # for struct_dir in sorted(glob(("zs/chai/struct_joint/PfTrpB*"))):
    #     dock_lib_parallel(
    #         struct_dir = struct_dir,
    #         dock_opt="substrate",
    #         score_only=False,
    #         rerun=False,
    #         cofactor_dets="cofactor",
    #         max_workers=32
    #     )
    
    # for struct_dir in sorted(glob(("zs/chai/struct_seperate/*"))):
    #     dock_lib_parallel(
    #         struct_dir = struct_dir,
    #         dock_opt="substrate",
    #         score_only=False,
    #         rerun=False,
    #         cofactor_dets="cofactor",
    #         max_workers=32
    #     )
    

    # for dock_opt in ["substrate", "all"]:
    #     for lib in sorted(glob(("data/meta/ParLQ-*.csv"))):
    #         if "TrpB" in lib:
    #             VinaLibDock(
    #                 input_csv = lib,
    #                 dock_opt=dock_opt,
    #             )
    #         else:
    #             VinaLibDock(
    #                 input_csv = lib,
    #                 dock_opt=dock_opt,
    #                 cofactor_dets="activated_carbene-cofactor",
    #                 redock= True,
    #                 max_workers=4
    #             )

    # for struct_dir in [
    #     "zs/chai/struct_joint/ParLQ",
    #     "zs/chai/struct_joint/Rma-CB",
    #     "zs/chai/struct_joint/Rma-CSi",
    # ]:
    #     dock_lib_parallel(
    #         struct_dir = struct_dir,
    #         dock_opt="substrate+carbene",
    #         score_only=True,
    #         rerun=False,
    #         cofactor_dets="cofactor",
    #         max_workers=8
    #     )
    
    # for struct_dir in [
    #     "zs/chai/struct_seperate/ParLQ",
    #     "zs/chai/struct_seperate/Rma-CB",
    #     "zs/chai/struct_seperate/Rma-CSi",
    # ]:
    #     dock_lib_parallel(
    #         struct_dir = struct_dir,
    #         dock_opt="substrate+carbene",
    #         score_only=True,
    #         rerun=False,
    #         cofactor_dets="cofactor",
    #         max_workers=8
    #     )
    

    # for struct_dir in [
    #     "zs/chai/struct_joint/ParLQ",
    #     "zs/chai/struct_joint/Rma-CB",
    #     "zs/chai/struct_joint/Rma-CSi",
    # ]:
    #     dock_lib_parallel(
    #         struct_dir = struct_dir,
    #         dock_opt="substrate+carbene",
    #         score_only=False,
    #         rerun=False,
    #         cofactor_dets="cofactor",
    #         max_workers=8
    #     )
    
    # for struct_dir in [
    #     "zs/chai/struct_seperate/ParLQ",
    #     "zs/chai/struct_seperate/Rma-CB",
    #     "zs/chai/struct_seperate/Rma-CSi",
    # ]:
    #     dock_lib_parallel(
    #         struct_dir = struct_dir,
    #         dock_opt="substrate+carbene",
    #         score_only=False,
    #         rerun=False,
    #         cofactor_dets="cofactor",
    #         max_workers=8
    #     )
    


    # for struct_dir in [
    #     "zs/af3/struct_joint/ParLQ",
    #     "zs/af3/struct_joint/Rma-CB",
    #     "zs/af3/struct_joint/Rma-CSi",
    # ]:
    #     dock_lib_parallel(
    #         struct_dir = struct_dir,
    #         dock_opt="substrate+carbene",
    #         score_only=True,
    #         rerun=False,
    #         cofactor_dets="cofactor",
    #         max_workers=8
    #     )
    
    # for struct_dir in [
    #     "zs/af3/struct_seperate/ParLQ",
    #     "zs/af3/struct_seperate/Rma-CB",
    #     "zs/af3/struct_seperate/Rma-CSi",
    # ]:
    #     dock_lib_parallel(
    #         struct_dir = struct_dir,
    #         dock_opt="substrate+carbene",
    #         score_only=True,
    #         rerun=False,
    #         cofactor_dets="cofactor",
    #         max_workers=8
    #     )
    

    # for struct_dir in [
    #     "zs/af3/struct_joint/ParLQ",
    #     "zs/af3/struct_joint/Rma-CB",
    #     "zs/af3/struct_joint/Rma-CSi",
    # ]:
    #     dock_lib_parallel(
    #         struct_dir = struct_dir,
    #         dock_opt="substrate+carbene",
    #         score_only=False,
    #         rerun=False,
    #         cofactor_dets="cofactor",
    #         max_workers=8
    #     )
    
    # for struct_dir in [
    #     "zs/af3/struct_seperate/ParLQ",
    #     "zs/af3/struct_seperate/Rma-CB",
    #     "zs/af3/struct_seperate/Rma-CSi",
    # ]:
    #     dock_lib_parallel(
    #         struct_dir = struct_dir,
    #         dock_opt="substrate+carbene",
    #         score_only=False,
    #         rerun=False,
    #         cofactor_dets="cofactor",
    #         max_workers=8
    #     )
    f.close()
    
    """
    VinaLibDock(ZSData):
    
    def __init__(
        self,
        input_csv: str,
        dock_opt: str,  #  ie "substrate", "joint", "all"
        cofactor_dets: str = "cofactor", # or inactivated_cofactor
        in_structure_dir: str = "data/structure",
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        fit_col_name: str = "fitness",
        output_dir: str = "zs/vina/apo",
        pH: float = 7.4,
        size_x: float = 20.0,
        size_y: float = 20.0,
        size_z: float = 20.0,
        num_modes: int = 9,
        exhaustiveness: int = 32,
        max_workers: int = 24,
        regen: bool = False,
        redock: bool = False
    )"""