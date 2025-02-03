"""Test the coves module."""

import sys
import os

from datetime import datetime

from REVIVAL.zs.coves import run_all_coves, append_all_coves_scores
from REVIVAL.util import checkNgen_folder

os.environ["CUDA_VISIBLE_DEVICES"] = "1"

if __name__ == "__main__":

    log_folder = checkNgen_folder("logs/zs/coves")

    # log outputs
    f = open(os.path.join(log_folder, f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out"), 'w')
    sys.stdout = f

    # run_all_coves(pattern="data/lib/*", withsub=False)
    # run_all_coves(pattern="data/lib/*", withsub=True)
    # run_all_coves(pattern="data/meta/*_scope.csv", withsub=False)
    # append_all_coves_scores(libs="data/meta/*",coves_dir="zs/coves/apo_clean/output/100")
    # append_all_coves_scores(libs="data/meta/*",coves_dir="zs/coves/sub_clean/output/100")

    f.close()

    """
    run_all_coves(
        pattern="data/lib/*",
        data_dir: str = "data",
        chain_number: str = "A",
        coves_dir="zs/coves",
        lmdb_dir: str = "lmdb",
        model_weight_path: str = "/disk2/fli/ddingding-CoVES/data/coves/res_weights/RES_1646945484.3030427_8.pt",
        dout: str = "zs/coves/output",
        n_ave: int = 100,
    )

    append_all_coves_scores(
        libs: list | str = "data/meta/scale2parent/*",
        input_dir: str = "data/meta/scale2parent",
        coves_dir: str = "zs/coves/output/100",
        t: float = 0.1,
    ) 
    """