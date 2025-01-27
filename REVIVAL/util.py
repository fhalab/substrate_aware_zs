"""
Util functions
"""

from __future__ import annotations

import re
import os
import json
import logging
import subprocess

import numpy as np

from Bio import SeqIO, pairwise2, PDB
from Bio.PDB import PDBParser, PDBIO, MMCIFParser

try:
    from rdkit import Chem
except:
    pass

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


def load_json(file_path: str) -> dict:
    """Load JSON content from a file."""
    with open(file_path, "r") as f:
        return json.load(f)


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


def get_protein_structure(file_path: str):

    """
    Get the parser based on the file extension

    Args:
    - file_path: str, the path to the file

    Returns:
    - PDBParser, the parser
    """

    # Check the file extension
    file_extension = os.path.splitext(file_path)[1].lower()

    if file_extension == ".pdb":
        parser = PDBParser(QUIET=True)  # QUIET suppresses warnings
    elif file_extension == ".cif":
        parser = MMCIFParser(QUIET=True)
    else:
        raise ValueError("Unsupported file format. Please use a .pdb or .cif file.")

    return parser.get_structure("protein", file_path)


def get_chain_ids(pdb_file_path: str) -> list:
    """
    Extract chain IDs from a given PDB file.

    Args:
    - pdb_file_path, str: The path to the PDB file.

    Returns:
    - list: A list of chain IDs in the PDB file.
    """
    # Parse the input PDB file

    structure = get_protein_structure(pdb_file_path)

    # Extract chain IDs
    chain_ids = set()
    for model in structure:
        for chain in model:
            chain_ids.add(chain.id)

    return sorted(list(chain_ids))


def get_chain_structure(input_file_path: str, output_file_path: str, chain_id: str):
    """
    Extract specified chains and save them to a new PDB file.

    Args:
        input_file_path (str): Path to the input PDB or CIF file.
        output_file_path (str): Path to save the output PDB file.
        chain_id (str): String of chain IDs to extract.
            ie. "A", "A,B", "A,B,C", etc.
    """
    # Check if input is a CIF file
    structure = get_protein_structure(input_file_path)

    # convert chain_id to a list
    chain_ids = chain_id.split(",")

    # Initialize PDBIO for writing the output
    io = PDBIO()
    extracted_chains = []

    if output_file_path:
        # Open the output file for writing
        with open(output_file_path, "w") as fout:
            for model in structure:
                for chain in model:
                    if chain.id in chain_ids:
                        # Save the chain to a temporary file-like object
                        io.set_structure(chain)
                        io.save(fout)  # Save chain to the output file
                        fout.write("\n")  # Add a newline for separation

        print(f"Chains {', '.join(chain_ids)} have been saved to {output_file_path}")
    else:
        for model in structure:
            for chain in model:
                if chain.id in chain_ids:
                    extracted_chains.append(chain)
        return extracted_chains


def calculate_chain_centroid(
    input_file: str, chain_ids
) -> np.ndarray:

    """
    Calculate the geometric center (centroid) of all atoms in the specified chain(s).

    Args:
        input_file (str): Path to the input PDB or CIF file.
        chain_ids (list of str): List of chain IDs to calculate the centroid for.

    Returns:
        tuple: The XYZ coordinates of the centroid.
    """

    # Parse the structure
    structure = get_protein_structure(input_file)

    coordinates = []
    chain_ids = [cid.upper() for cid in chain_ids]  # Ensure chain IDs are uppercase

    for model in structure:
        for chain in model:
            if chain.id.upper() in chain_ids:
                for residue in chain:
                    for atom in residue:
                        coordinates.append(atom.coord)

    # Calculate centroid
    if coordinates:
        centroid = np.mean(coordinates, axis=0)
        return np.array(centroid).flatten()
    else:
        raise ValueError(f"No atoms found for the specified chain(s): {chain_ids}")


def calculate_ligand_centroid(pdb_file: str, ligand_info: list):
    """
    Calculates the centroid of a given list of atoms specified by chain, residue, and atom names.

    Args:
        pdb_file (str): Path to the PDB file.
        ligand_info (list of tuples): List of atoms specified as (chain_id, residue_name, atom_name).

    Returns:
        tuple: Centroid coordinates as (x, y, z).
    """

    structure = get_protein_structure(pdb_file)

    atom_coords = []

    for chain_id, residue_name, atom_name in ligand_info:
        for model in structure:
            try:
                chain = model[chain_id]
                for residue in chain:
                    # Match residue name (flexible: partial match, case insensitive)
                    if residue_name.lower() in residue.resname.lower():
                        # Match atom name (flexible: ignore underscores and case differences)
                        for atom in residue:
                            if atom_name.replace("_", "").lower() == atom.name.replace("_", "").lower():
                                atom_coords.append(atom.coord)
            except KeyError:
                print(f"Chain {chain_id} not found in the structure.")
                continue

    if not atom_coords:
        raise ValueError("No matching atoms found in the structure.")

    # Calculate centroid
    return np.mean(atom_coords, axis=0)
    

def modify_PDB_chain(
    input_file_path: str,
    output_file_path: str,
    original_chain_id: str,
    modified_chain_id: str,
):
    """
    Modify chain ID in a PDB file and save it to a new file.

    Parameters:
    input_file_path (str): Path to the input PDB file.
    output_file_path (str): Path to the output PDB file.
    original_chain_id (str): The chain ID to be replaced.
    modified_chain_id (str): The new chain ID.
    """

    # Parse the input PDB file
    structure = get_protein_structure(input_file_path)

    # Iterate over the chains and replace the original chain ID with the modified chain ID
    for model in structure:
        for chain in model:
            if chain.id == original_chain_id:
                chain.id = modified_chain_id

    # Save the modified structure to the output PDB file
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file_path)

    print(
        f"Chain {original_chain_id} has been replaced with {modified_chain_id} in {output_file_path}"
    )


def cif2pdbwobabel(cif_file: str, pdb_file: str):

    cmd = f"obabel {cif_file} -O {pdb_file} --remove HOH"
    subprocess.run(cmd, shell=True)
    

def convert_cif_to_pdb(cif_file: str, pdb_file: str = "", ifsave: bool = True):
    """
    Converts a CIF file to PDB format while preserving as much structural data as possible.

    Args:
    - cif_file (str): Path to the input CIF file.
    - pdb_file (str): Path to the output PDB file.
    """

    # Parse the CIF file
    structure = get_protein_structure(cif_file)

    # Create a PDBIO object to write the structure to PDB format
    io = PDBIO()

    # Set the structure to write
    io.set_structure(structure)

    if ifsave:
        # Save the structure to the PDB file
        io.save(pdb_file)

        print(f"Converted {cif_file} to {pdb_file}")
    else:
        return structure


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
    structure = get_protein_structure(input_pdb)

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


def replace_residue_names_auto(
    input_file: str, 
    output_file: str,
    residue_prefix: str = "LIG",
    new_residue: str = "LIG"):
    """
    Automatically detect and replace residue names in a PDB file that match a specific prefix.

    Args:
        input_file (str): Path to the input PDB file.
        output_file (str): Path to save the modified PDB file.
        residue_prefix (str): Prefix of residue names to replace (e.g., "LIG").
        new_residue (str): New residue name to replace with.
    """

    # Convert CIF to PDB if necessary
    if input_file.lower().endswith(".cif"):
        pdb_file = os.path.splitext(input_file)[0] + ".pdb"
        print(f"Converting CIF to PDB: {input_file} -> {pdb_file}")
        
        # Parse CIF and write as PDB
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("structure", input_file)
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_file)
        input_file = pdb_file  # Use the converted PDB file as input

    detected_residues = set()  # To store dynamically detected residue names
    pattern = re.compile(f"^{residue_prefix}_\\w$")  # Regex to detect residues like LIG_B, LIG_C

    with open(input_file, "r") as infile:
        lines = infile.readlines()

    # First pass: Detect residue names dynamically
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            long_res_name = line[17:22]
            if pattern.match(long_res_name):
                detected_residues.add(long_res_name)

    if len(detected_residues) > 0:
        print(f"Detected residues to replace: {detected_residues}")

    # batch replace detected residues with new residue name
    with open(output_file, "w") as outfile:
        for line in lines:
            if line.startswith(("ATOM", "HETATM")):
                long_res_name = line[17:22]
                if long_res_name in detected_residues:
                    line = line.replace(long_res_name, new_residue)

                # Clean up atom names by removing underscores
                atom_name = line[12:16].strip()
                cleaned_atom_name = atom_name.replace("_", "")
                line = line[:12] + cleaned_atom_name.ljust(4) + line[16:]

            outfile.write(line)


from Bio.PDB import PDBParser, PDBIO, Select

class ResidueRemover(Select):
    def __init__(self, residues_to_remove):
        """
        residues_to_remove: list of residue names to remove (e.g., ['PLS', 'NA'])
        """
        self.residues_to_remove = residues_to_remove

    def accept_residue(self, residue):
        """Return False if the residue should be removed."""
        if residue.get_resname() in self.residues_to_remove:
            return False
        return True


def remove_residues_from_pdb(
    input_pdb: str, 
    output_pdb: str, 
    residues_to_remove: list
):
    """
    Remove specified residues from a PDB file.

    Args:
        input_pdb (str): Path to the input PDB file.
        output_pdb (str): Path to save the modified PDB file.
        residues_to_remove (list): List of residue names to remove (e.g., ['PLS', 'NA']).
    """
    # Parse the input PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb)

    # Save the modified structure to a new PDB file
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, select=ResidueRemover(residues_to_remove))


def remove_hetatm(input_pdb: str, output_pdb: str):
    """
    Remove all HETATM entries from a PDB file.

    Args:
        input_pdb (str): Path to the input PDB file.
        output_pdb (str): Path to save the modified PDB file.
    """
    with open(input_pdb, "r") as infile, open(output_pdb, "w") as outfile:
        for line in infile:
            # Write the line to the output file if it doesn't start with "HETATM"
            if not line.startswith("HETATM"):
                outfile.write(line)


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


def er2ee(er: str) -> float:

    """
    Convert enantiomeric ratio to enantiomeric excess
    > 99.9:0.1 is considered as 99.8% ee
    """

    pdt1, pdt2 = map(
        float, er.replace(">", "").split(":")
    )  # Split the ratio into major and minor components
    return (pdt1 - pdt2) / (pdt1 + pdt2) * 100  # Apply the EE formula


def canonicalize_smiles(smiles_string: str) -> str:

    """
    A function to canonicalize a SMILES string.

    Args:
    - smiles_string (str): The input SMILES string.

    Returns:
    - str: The canonicalized SMILES string.
    """

    molecule = Chem.MolFromSmiles(smiles_string)
    if molecule:
        canonical_smiles = Chem.MolToSmiles(molecule, canonical=True)
        return canonical_smiles


def smiles2mol(smiles_string: str):
    """
    A function to convert a SMILES string to an RDKit molecule object.

    Args:
    - smiles_string (str): The input SMILES string.

    Returns:
    - RDKit molecule object.
    """

    mol = Chem.MolFromSmiles(smiles_string)
    if mol:
        mol = Chem.AddHs(mol)
        return mol
    else:
        raise ValueError("Invalid SMILES string.")


def protonate_smiles(smiles: str, pH: float) -> str:
    """
    Protonate SMILES string with OpenBabel at given pH

    :param smiles: SMILES string of molecule to be protonated
    :param pH: pH at which the molecule should be protonated
    :return: SMILES string of protonated structure
    """

    # cmd list format raises errors, therefore one string
    cmd = f'obabel -:"{smiles}" -ismi -ocan -p{pH}'
    cmd_return = subprocess.run(cmd, capture_output=True, shell=True)
    output = cmd_return.stdout.decode("utf-8")
    logging.debug(output)

    if cmd_return.returncode != 0:
        print("WARNING! COULD NOT PROTONATE")
        return None

    return output.strip()


def protonate_oxygen(smiles: str) -> str:
    """
    Protonate all [O-] groups in a SMILES string.

    :param smiles: Input SMILES string with [O-] groups.
    :return: Protonated SMILES string with [OH] instead of [O-].
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    # Add hydrogens explicitly
    mol = Chem.AddHs(mol)

    # Iterate over atoms to find [O-] and adjust charges
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "O" and atom.GetFormalCharge() == -1:
            # Set the charge to neutral
            atom.SetFormalCharge(0)
            # Adjust the number of implicit hydrogens
            atom.SetNumExplicitHs(1)

    # Update the molecule
    Chem.SanitizeMol(mol)

    # Generate the protonated SMILES
    protonated_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    return protonated_smiles


def add_hydrogens_to_smiles(smiles: str) -> str:
    """
    Add explicit hydrogens to a molecule represented by a SMILES string.

    Args:
        smiles (str): Input SMILES string.

    Returns:
        str: SMILES string with explicit hydrogens.
    """
    mol = Chem.MolFromSmiles(smiles)  # Parse SMILES to RDKit molecule
    if not mol:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    
    mol_with_h = Chem.AddHs(mol)  # Add explicit hydrogens
    smiles_with_h = Chem.MolToSmiles(mol_with_h)  # Convert back to SMILES
    
    return smiles_with_h


def run_sh_command(command, working_dir=None):
    """
    Run a shell command in a specified directory.

    Parameters:
    - command: The shell command to execute.
    - working_dir: The directory where the command will be run (optional).
    """
    try:
        subprocess.run(command, shell=True, check=True, cwd=working_dir)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")