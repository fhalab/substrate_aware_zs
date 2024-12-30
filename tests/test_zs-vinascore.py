"""
Test the vina docking and extracting score.
"""

import sys
import os

from glob import glob

from datetime import datetime

from REVIVAL.zs.vina import run_parse_vina_results
from REVIVAL.util import checkNgen_folder

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

    run_parse_vina_results(
        pattern="data/meta/not_scaled/PfTrpB-*.csv",
        cofactor_type="inactivated-cofactor",
    )

    # def run_parse_vina_results(
    #     pattern: Union[str, list] = "data/meta/not_scaled/*.csv",
    #     cofactor_type: str = "",
    #     kwargs: dict = {},
    # )

    f.close()