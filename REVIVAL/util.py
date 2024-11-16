"""
Util functions
"""

from __future__ import annotations

import os

from Bio import SeqIO, pairwise2, PDB
from Bio.PDB import PDBParser, PDBIO, MMCIFParser
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


def get_chain(input_file_path: str, output_file_path: str, chain_id: str):

    """
    Get the chain given ID
    """

    # Parse the input PDB file
    parser = PDBParser()
    structure = parser.get_structure("protein", input_file_path)

    # Iterate over the chains and replace the original chain ID with the modified chain ID
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                io = PDBIO()
                io.set_structure(chain)
                io.save(output_file_path)

    print(f"Chain {chain_id} has been saved to {output_file_path}")


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


def convert_cif_to_pdb(cif_file: str, pdb_file: str):
    """
    Converts a CIF file to PDB format while preserving as much structural data as possible.

    Args:
    - cif_file (str): Path to the input CIF file.
    - pdb_file (str): Path to the output PDB file.
    """
    # Create a CIF parser object
    parser = MMCIFParser(QUIET=True)
    
    # Parse the CIF file
    structure = parser.get_structure('structure', cif_file)

    # Create a PDBIO object to write the structure to PDB format
    io = PDBIO()

    # Set the structure to write
    io.set_structure(structure)

    # Save the structure to the PDB file
    io.save(pdb_file)

    print(f"Converted {cif_file} to {pdb_file}")


def chop_pdb(
    input_pdb: str, output_pdb: str, start_resid: int, end_resid: int, chain_id: str
) -> None:
    """
    A function for chopping a pdb file to a specific chain and residue range

    Args:
    - input_pdb: str, path to the input pdb file
    - output_pdb: str, path to the output pdb file
    - start_resid: int, starting residue ID
    - end_resid: int, ending residue ID
    - chain_id: str, chain ID
    """

    # Initialize the parser and structure
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("structure", input_pdb)

    # Initialize the writer
    io = PDB.PDBIO()

    # Define a select class to filter the residues in the specific chain
    class ResidueSelect(PDB.Select):
        def accept_residue(self, residue):
            # Only accept residues in the specified chain with a residue ID greater than or equal to start_resid
            if (
                residue.parent.id == chain_id
                and residue.id[1] >= start_resid
                and residue.id[1] <= end_resid
            ):
                return True
            return False

    # Save the chopped structure to the output file
    io.set_structure(structure)
    io.save(output_pdb, ResidueSelect())

    print(
        f"Saved chopped structure starting from residue {start_resid} in chain {chain_id} to {output_pdb}"
    )


def pdb2seq(pdb_file_path: str, chain_id: str = "A") -> str:

    """
    A function for extracting chain in string format from pdb

    Args:
    - pdb_file_path: str,
    - chain_id: str = "A"
    """

    chains = {
        record.id: record.seq for record in SeqIO.parse(pdb_file_path, "pdb-atom")
    }

    return str(chains[[chain for chain in chains.keys() if chain_id in chain][0]])


def find_missing_str(longer: str, shorter: str) -> [str, str]:
    """
    A function for finding the missing part of a string

    Args:
    - longer: str, longer string
    - shorter: str, shorter string

    Returns:
    - part_before: str, part of the longer string before the shorter
    - part_after: str, part of the longer string after the shorter
    """
    # Find the start index of the shorter in the longer string
    start_index = longer.find(shorter)

    # If the shorter is not found, return the longer string as the "missing" part
    if start_index == -1:
        return "", ""

    # Find the end index of the shorter
    end_index = start_index + len(shorter)

    # Extract parts of the longer string that are not the shorter
    part_before = longer[:start_index]
    part_after = longer[end_index:]

    return part_before, part_after


def alignmutseq2pdbseq(mut_seq: str, pdb_seq: str) -> list[int]:
    """
    A function for aligning mutation sequence to pdb sequence and
    return the indices of the aligned sequence so that the mutation
    sequence can be trimmed to the lenght of the pdb sequence

    Args:
    - mut_seq: str, mutation sequence
    - pdb_seq: str, pdb sequence

    Returns:
    - list[int], start and end indices of the aligned sequence
    - pdb_seq: str, aligned pdb sequence
    """

    # Define a custom scoring function so that X is aligned with anything
    def custom_match_function(x, y):
        if x == "X" or y == "X":
            return 2  # High score for aligning X with anything
        elif x == y:
            return 2  # Match score
        else:
            return -1  # Mismatch score

    _, aligned_pdb_seq, _, _, _ = pairwise2.align.globalcs(
        mut_seq, pdb_seq, custom_match_function, -0.5, -0.1
    )[0]

    return [
        aligned_pdb_seq.find(aligned_pdb_seq.replace("-", "")[:1]),
        aligned_pdb_seq.rfind(aligned_pdb_seq.replace("-", "")[-1]),
    ], aligned_pdb_seq
