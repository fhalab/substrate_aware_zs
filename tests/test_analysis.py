"""Test the preprocess module."""

import sys
import os

from datetime import datetime

from REVIVAL.analysis import train_test_all
from REVIVAL.util import checkNgen_folder


if __name__ == "__main__":

    log_folder = checkNgen_folder("logs/lincomb")

    # log outputs
    f = open(os.path.join(log_folder, f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out"), 'w')
    sys.stdout = f


    train_test_all(
        pattern="/disk2/fli/REVIVAL2/zs/comb/*.csv", output_dir="zs/lincomb"
    )

    f.close()

    """
    train_test_all(
        pattern="/disk2/fli/REVIVAL2/zs/comb/*.csv", output_dir="zs/lincomb"
    )
    """


