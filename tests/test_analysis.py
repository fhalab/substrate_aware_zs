"""Test the preprocess module."""

import sys
import os

from datetime import datetime

from substrate_aware.analysis import train_test_all, plot_all_metrics, get_minimal_comb # process_and_save_metrics
from substrate_aware.util import checkNgen_folder


if __name__ == "__main__":

    log_folder = checkNgen_folder("logs/analysis")

    # log outputs
    f = open(os.path.join(log_folder, f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out"), 'w')
    sys.stdout = f

    # get_minimal_comb(pattern="zs/comb/*.csv")

    # train_test_all(
    #     pattern="/disk2/fli/substrate_aware/zs/comb/minimal/*.csv", output_dir="zs/lincomb"
    # )
    # train_test_all(
    #     pattern="/disk2/fli/substrate_aware/zs/comb/minimal/*.csv",
    #     output_dir="zs/lincomb",
    #     sele_cols=["EVmutation", "AF3"]
    # )

    # train_test_all(
    #     pattern="/disk2/fli/substrate_aware/zs/comb/minimal/*.csv",
    #     output_dir="zs/lincomb",
    #     sele_cols=["EVmutation", "AF3", "ESM-IF"]
    # )

    # train_test_all(
    #     pattern="/disk2/fli/substrate_aware/zs/comb/minimal/*.csv",
    #     output_dir="zs/lincomb",
    #     sele_cols=["EVmutation", "ESM-IF"]
    # )
    # train_test_all(
    #     pattern="/disk2/fli/substrate_aware/zs/comb/minimal/*.csv",
    #     output_dir="zs/lincomb",
    #     sele_cols=["AF3", "ESM-IF"]
    # )


    # process_and_save_metrics(
    #     input_dir="zs/comb/minimal", 
    #     output_dir="zs/metrics"
    #     )

    # plot_all_metrics(
    #     input_dir="zs/metrics",
    #     output_dir="figs/metrics"
    #     )

    f.close()

