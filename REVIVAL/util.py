"""
Util functions
"""

from __future__ import annotations

import os

from Bio import SeqIO
# import pickle

# import numpy as np
# import pandas as pd

# from sklearn.metrics import ndcg_score

def checkNgen_folder(folder_path: str) -> str:

    """
    Check if the folder and its subfolder exists
    create a new directory if not
    Args:
    - folder_path: str, the folder path
    """

    # if input path is file
    if bool(os.path.splitext(folder_path)[1]):
        folder_path = os.path.dirname(folder_path)

    split_list = os.path.normpath(folder_path).split("/")
    for p, _ in enumerate(split_list):
        subfolder_path = "/".join(split_list[: p + 1])
        if not os.path.exists(subfolder_path):
            print(f"Making {subfolder_path} ...")
            os.mkdir(subfolder_path)
    return folder_path


def get_file_name(file_path: str) -> str:

    """
    Extract file name without the extension
    Args:
    - file_path: str, ie. data/graph_nx/Tm9D8s/Tm9D8s_3siteA_fixed/WT.pdb
    Returns:
    - str, ie WT
    """

    return os.path.splitext(os.path.basename(file_path))[0]


def get_dir_name(file_path: str) -> str:

    """
    Extract dir name
    Args:
    - file_path: str, ie. data/graph_nx/Tm9D8s/Tm9D8s_3siteA_fixed/WT.pdb
    Returns:
    - str, ie Tm9D8s_3siteA_fixed
    """

    return os.path.basename(os.path.dirname(file_path))


def read_parent_fasta(fasta_path: str) -> str:

    """
    Read the parent fasta file
    
    Args:
    - fasta_path: str, the path to the fasta file

    Returns:
    - str, the parent sequence
    """

    # Parse the fasta file
    sequences = list(SeqIO.parse(fasta_path, "fasta"))
    
    # Assert that there is exactly one sequence in the file
    assert len(sequences) == 1, "FASTA file contains more than one sequence."

    # Return the sequence as a string
    return str(sequences[0].seq)