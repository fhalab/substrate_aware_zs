"""
Util functions
"""

from __future__ import annotations

import os

from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO
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



def get_chain_ids(pdb_file_path: str) -> list:
    """
    Extract chain IDs from a given PDB file.

    Args:
    - pdb_file_path (str): The path to the PDB file.

    Returns:
    - list: A list of chain IDs in the PDB file.
    """
    parser = PDBParser(QUIET=True)  # QUIET=True suppresses warnings
    structure = parser.get_structure("structure", pdb_file_path)

    # Extract chain IDs
    chain_ids = set()
    for model in structure:
        for chain in model:
            chain_ids.add(chain.id)

    return list(chain_ids)


def modify_PDB_chain(input_file_path: str, output_file_path: str, original_chain_id: str, modified_chain_id: str):
    """
    Modify chain ID in a PDB file and save it to a new file.

    Parameters:
    input_file_path (str): Path to the input PDB file.
    output_file_path (str): Path to the output PDB file.
    original_chain_id (str): The chain ID to be replaced.
    modified_chain_id (str): The new chain ID.
    """

    # Parse the input PDB file
    parser = PDBParser()
    structure = parser.get_structure("protein", input_file_path)

    # Iterate over the chains and replace the original chain ID with the modified chain ID
    for model in structure:
        for chain in model:
            if chain.id == original_chain_id:
                chain.id = modified_chain_id

    # Save the modified structure to the output PDB file
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file_path)

    print(f"Chain {original_chain_id} has been replaced with {modified_chain_id} in {output_file_path}")
