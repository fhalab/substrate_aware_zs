"""
Test the vina docking and extracting score.
"""

import sys
import os

from glob import glob

from datetime import datetime

from substrate_aware.zs.vina import run_parse_vina_results
from substrate_aware.util import checkNgen_folder

if __name__ == "__main__":

    os.environ["CUDA_VISIBLE_DEVICES"] = "1"

    # log outputs
    f = open(
        os.path.join(
            checkNgen_folder("logs/zs/vina/dock"),
            f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out",
        ),
        "w",
    )
    sys.stdout = f

    # run_parse_vina_results(
    #     pattern="data/meta/ParLQ-*.csv", # "data/meta/*.csv",
    #     dock_opt="substrate",
    #     score_only=False,
    #     vina_struct_dir="zs/vina/apo/struct",
    # )


    # run_parse_vina_results(
    #     pattern="data/meta/Rma-CB.csv", # "data/meta/*.csv",
    #     dock_opt="substrate",
    #     score_only=False,
    #     vina_struct_dir="zs/vina/apo/struct",
    # )

    # run_parse_vina_results(
    #     pattern="data/meta/Rma-CB.csv", # "data/meta/*.csv",
    #     dock_opt="all",
    #     score_only=False,
    #     vina_struct_dir="zs/vina/apo/struct",
    # )

    # run_parse_vina_results(
    #     pattern="data/meta/*.csv",
    #     dock_opt="substrate",
    #     score_only=True,
    #     vina_struct_dir="zs/vina/chai/struct_joint",
    # )

    # run_parse_vina_results(
    #     pattern="data/meta/ParLQ-*.csv",
    #     dock_opt="substrate",
    #     score_only=True,
    #     vina_struct_dir="zs/vina/chai/struct_seperate",
    # )

    # run_parse_vina_results(
    #     pattern="data/meta/*.csv",
    #     dock_opt="substrate",
    #     score_only=True,
    #     vina_struct_dir="zs/vina/af3/struct_joint",
    # )

    # run_parse_vina_results(
    #     pattern="data/meta/ParLQ-*.csv",
    #     dock_opt="substrate",
    #     score_only=True,
    #     vina_struct_dir="zs/vina/af3/struct_seperate",
    # )


    # run_parse_vina_results(
    #     pattern="data/meta/*.csv",
    #     dock_opt="substrate",
    #     score_only=False,
    #     vina_struct_dir="zs/vina/chai/struct_joint",
    # )

    # run_parse_vina_results(
    #     pattern="data/meta/*.csv",
    #     dock_opt="substrate",
    #     score_only=False,
    #     vina_struct_dir="zs/vina/chai/struct_seperate",
    # )

    # run_parse_vina_results(
    #     pattern="data/meta/*.csv",
    #     dock_opt="substrate",
    #     score_only=False,
    #     vina_struct_dir="zs/vina/af3/struct_joint",
    # )

    # run_parse_vina_results(
    #     pattern="data/meta/*.csv",
    #     dock_opt="substrate",
    #     score_only=False,
    #     vina_struct_dir="zs/vina/af3/struct_seperate",
    # )

    # run_parse_vina_results(
    #     pattern=[
    #         "data/meta/ParLQ.csv",
    #         "data/meta/Rma-CB.csv",
    #         "data/meta/Rma-CSi.csv",
    #     ],
    #     dock_opt="substrate+carbene",
    #     score_only=False,
    #     vina_struct_dir="zs/vina/chai/struct_seperate",
    # )

    # run_parse_vina_results(
    #     pattern=[
    #         "data/meta/ParLQ.csv",
    #         "data/meta/Rma-CB.csv",
    #         "data/meta/Rma-CSi.csv",
    #     ],
    #     dock_opt="substrate+carbene",
    #     score_only=False,
    #     vina_struct_dir="zs/vina/chai/struct_joint",
    # )
    
    # run_parse_vina_results(
    #     pattern=[
    #         "data/meta/ParLQ.csv",
    #         "data/meta/Rma-CB.csv",
    #         "data/meta/Rma-CSi.csv",
    #     ],
    #     dock_opt="substrate+carbene",
    #     score_only=True,
    #     vina_struct_dir="zs/vina/chai/struct_joint",
    # )
    

    # run_parse_vina_results(
    #     pattern=[
    #         "data/meta/ParLQ.csv",
    #         "data/meta/Rma-CB.csv",
    #         "data/meta/Rma-CSi.csv",
    #     ],
    #     dock_opt="substrate+carbene",
    #     score_only=True,
    #     vina_struct_dir="zs/vina/chai/struct_seperate",
    # )
    
    
    # run_parse_vina_results(
    #     pattern=[
    #         "data/meta/ParLQ.csv",
    #         "data/meta/Rma-CB.csv",
    #         "data/meta/Rma-CSi.csv",
    #     ],
    #     dock_opt="substrate+carbene",
    #     score_only=True,
    #     vina_struct_dir="zs/vina/af3/struct_joint",
    # )

    # run_parse_vina_results(
    #     pattern=[
    #         "data/meta/ParLQ.csv",
    #         "data/meta/Rma-CB.csv",
    #         "data/meta/Rma-CSi.csv",
    #     ],
    #     dock_opt="substrate+carbene",
    #     score_only=False,
    #     vina_struct_dir="zs/vina/af3/struct_joint",
    # )

    # run_parse_vina_results(
    #     pattern=[
    #         "data/meta/ParLQ.csv",
    #         "data/meta/Rma-CB.csv",
    #         "data/meta/Rma-CSi.csv",
    #     ],
    #     dock_opt="substrate+carbene",
    #     score_only=True,
    #     vina_struct_dir="zs/vina/af3/struct_seperate",
    # )

    # run_parse_vina_results(
    #     pattern=[
    #         "data/meta/ParLQ.csv",
    #         "data/meta/Rma-CB.csv",
    #         "data/meta/Rma-CSi.csv",
    #     ],
    #     dock_opt="substrate+carbene",
    #     score_only=False,
    #     vina_struct_dir="zs/vina/af3/struct_seperate",
    # )

    # run_parse_vina_results(
    # dock_opt: str,  #  ie "substrate",
    # score_only: bool,  # = True,
    # vina_struct_dir: str, # = "zs/vina/chai/struct_joint",
    # pattern: Union[str, list] = "data/meta/*.csv",
    # kwargs: dict = {},

    f.close()