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

    for chai_dir in [
        # "zs/chai/mut_structure/PfTrpB-4bromo_cofactor",
        # "zs/chai/mut_structure/PfTrpB-4cyano_cofactor",
        # "zs/chai/mut_structure/PfTrpB-5bromo_cofactor",
        "zs/chai/mut_structure/PfTrpB-7iodo_cofactor",
        "zs/chai/mut_structure/PfTrpB-7bromo_cofactor",
        "zs/chai/mut_structure/PfTrpB-6chloro_cofactor",
        "zs/chai/mut_structure/PfTrpB-5iodo_cofactor",
        "zs/chai/mut_structure/PfTrpB-5cyano_cofactor",
        "zs/chai/mut_structure/PfTrpB-5chloro_cofactor",
        
    ]:
        dock_lib_parallel(
            chai_dir=chai_dir, 
            cofactor_type="cofactor",
            freeze_opt="cofactor",
            max_workers=16
            )

    # for chai_dir in sorted(glob("zs/chai/mut_structure/*_inactivated-cofactor")):
    #     dock_lib_parallel(chai_dir, "inactivated-cofactor")

    # dock_lib_parallel(
    #     chai_dir: str,
    #     cofactor_type: str,
    #     residues: list = None,
    #     substrate_chain_ids: Union[list, str] = "B",
    #     freeze_opt: str = None,
    #     vina_dir: str = "vina",
    #     pH: float = 7.4,
    #     method="vina",
    #     size_x=10.0,
    #     size_y=10.0,
    #     size_z=10.0,
    #     num_modes=9,
    #     exhaustiveness=32,
    #     regen=False,
    #     rerun=False,
    #     max_workers=24,  # Number of parallel workers
    # ):

    f.close()