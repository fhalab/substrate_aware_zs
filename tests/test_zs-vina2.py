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
        "zs/chai/mut_structure/PfTrpB-4bromo_inactivated-cofactor", 
        "zs/chai/mut_structure/PfTrpB-4cyano_inactivated-cofactor",
        "zs/chai/mut_structure/PfTrpB-5chloro_inactivated-cofactor",
        "zs/chai/mut_structure/PfTrpB-5bromo_inactivated-cofactor", 
        "zs/chai/mut_structure/PfTrpB-6chloro_inactivated-cofactor",
        "zs/chai/mut_structure/PfTrpB-7bromo_inactivated-cofactor",
        "zs/chai/mut_structure/PfTrpB-7iodo_inactivated-cofactor",
        "zs/chai/mut_structure/PfTrpB-7methyl_inactivated-cofactor",
        "zs/chai/mut_structure/PfTrpB-56chloro_inactivated-cofactor",
        "zs/chai/mut_structure/PfTrpB-5cyano_inactivated-cofactor",
    ]:
        dock_lib_parallel(chai_dir, "inactivated-cofactor")


    # def dock_lib(
    #     chai_dir: str,
    #     cofactor_type: str,
    #     vina_dir: str = "vina",
    #     pH: float = 7.4,
    #     method='vina',
    #     size_x=15.0,
    #         size_y=15.0, 
    #         size_z=15.0,
    #         num_modes=9, # Dunno check vina docks using the defaut
    #         exhaustiveness=32 
    # ):

    f.close()