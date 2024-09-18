"""Test the preprocess module."""

import sys
import os

from datetime import datetime

from REVIVAL.preprocess import preprocess_all
from REVIVAL.util import checkNgen_folder


if __name__ == "__main__":

    log_folder = checkNgen_folder("logs/preprocess")

    # log outputs
    f = open(os.path.join(log_folder, f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out"), 'w')
    sys.stdout = f

    preprocess_all(scale_fit="parent")

    f.close()

    """
    preprocess_all(
        input_pattern: str = "data/lib/*.csv",
        scale_fit: str = "parent",
        combo_col_name: str = "AAs",
        mut_col_name: str = "muts",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        seq_dir: str = "data/seq",
        output_dir: str = "data/meta",
    )
    """