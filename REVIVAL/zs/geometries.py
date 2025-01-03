"""
A script for calcalating the bond distances from generated structures
"""

from __future__ import annotations

import os
import re
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, PDBIO

from REVIVAL.preprocess import ZSData
from REVIVAL.util import checkNgen_folder


def get_atom_with_variations(residue, atom_name):
    """
    Attempt to retrieve an atom from a residue, trying multiple variations of the atom name.

    Args:
        residue: A Bio.PDB.Residue object.
        atom_name (str): The atom name to search for.

    Returns:
        Atom object if found.

    Raises:
        KeyError: If no matching atom name is found.
    """

    variations = [
        atom_name,
        atom_name.replace("_", ""),  # Remove underscores
        f"{atom_name[0]}_{atom_name[1:]}",  # Add underscore after the first character
    ]
    for variation in variations:
        try:
            return residue[variation]
        except KeyError:
            continue
    raise KeyError(f"Atom {atom_name} or its variations not found in residue {residue}.")


# Match residues dynamically
def find_residue_by_id(chain, target_res_id):
    """
    Find a residue in a chain based on its ID.

    Args:
        chain: Biopython Chain object.
        target_res_id: Residue ID (sequence number) to search for.

    Returns:
        Residue object if found.
    
    Raises:
        ValueError: If the residue is not found in the chain.
    """

    for residue in chain.get_residues():
        if residue.id[1] == target_res_id:
            return residue
    raise ValueError(f"Residue with ID {target_res_id} not found in chain.")

    
def replace_residue_names_auto(input_file, output_file, residue_prefix="LIG", new_residue="LIG"):
    """
    Automatically detect and replace residue names in a PDB file that match a specific prefix.

    Args:
        input_file (str): Path to the input PDB file.
        output_file (str): Path to save the modified PDB file.
        residue_prefix (str): Prefix of residue names to replace (e.g., "LIG").
        new_residue (str): New residue name to replace with.
    """

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

    print(f"Detected residues to replace: {detected_residues}")

    # batch replace detected residues with new residue name
    with open(output_file, "w") as outfile:
        for line in lines:
            if line.startswith(("ATOM", "HETATM")):
                long_res_name = line[17:22]
                if long_res_name in detected_residues:
                    line = line.replace(long_res_name, new_residue)
                outfile.write(line)


def get_covalent_neighbors(atom, residue, bond_distance_cutoff=1.6):
    """
    Get covalently bonded neighbors of an atom in a residue based on distance.

    Args:
        atom: Biopython Atom object for the input atom.
        residue: Biopython Residue object containing the atom.
        bond_distance_cutoff (float): Maximum distance for covalent bonds in Ångströms.

    Returns:
        List: Neighboring atoms covalently bonded to the input atom.
    """

    atom_coord = np.array(atom.coord)
    neighbors = []
    for neighbor in residue.get_atoms():
        if neighbor != atom:  # Exclude the atom itself
            neighbor_coord = np.array(neighbor.coord)
            distance = np.linalg.norm(atom_coord - neighbor_coord)
            if distance <= bond_distance_cutoff:
                neighbors.append(neighbor)
    return neighbors


def calculate_hydrogen_position(atom, neighbors, bond_length=1.0):
    """
    Calculate the coordinates of a hydrogen atom based on the input atom's geometry.

    Args:
        atom: Biopython Atom object for the input atom.
        neighbors: List of neighboring Biopython Atom objects.
        bond_length (float): Bond length for the hydrogen atom in Ångströms.

    Returns:
        np.ndarray: Hydrogen atom coordinates.
    """

    atom_coord = np.array(atom.coord)

    if len(neighbors) == 1:  # SP hybridization (linear)
        neighbor_coord = np.array(neighbors[0].coord)
        direction = atom_coord - neighbor_coord
        direction /= np.linalg.norm(direction)  # Normalize
        hydrogen_coord = atom_coord + direction * bond_length

    elif len(neighbors) == 2:  # SP2 hybridization (planar)
        neighbor_coords = [np.array(neighbor.coord) for neighbor in neighbors]
        v1 = neighbor_coords[0] - atom_coord
        v2 = neighbor_coords[1] - atom_coord
        v1 /= np.linalg.norm(v1)
        v2 /= np.linalg.norm(v2)
        in_plane_direction = -(v1 + v2)  # Opposite direction within the plane
        in_plane_direction /= np.linalg.norm(in_plane_direction)
        hydrogen_coord = atom_coord + in_plane_direction * bond_length

    elif len(neighbors) == 3:  # SP3 hybridization (tetrahedral)
        # TODO test this case
        neighbor_coords = [np.array(neighbor.coord) for neighbor in neighbors]
        centroid = np.mean(neighbor_coords, axis=0)
        direction = atom_coord - centroid
        direction /= np.linalg.norm(direction)  # Normalize
        hydrogen_coord = atom_coord + direction * bond_length

    else:
        raise ValueError(f"Unsupported geometry for atom with {len(neighbors)} neighbors.")

    print(f"Calculated hydrogen position from {atom} based on {neighbors}: {hydrogen_coord}")

    return hydrogen_coord


def measure_bond_distance(
    structure_file, 
    chain_id_1, 
    res_id_1, 
    atom_name_1, 
    chain_id_2, 
    res_id_2, 
    atom_name_2, 
    add_hydrogen_to_1=False, 
    add_hydrogen_to_2=False):
    """
    Measure the bond distance between two atoms or hydrogens attached to them in a PDB or CIF file.

    Args:
        structure_file (str): Path to the PDB or CIF file.
        chain_id_1 (str): Chain ID where the first atom is located.
        res_id_1 (tuple): Tuple of (residue sequence number, insertion code) for the first atom.
        atom_name_1 (str): Name of the first atom.
        chain_id_2 (str): Chain ID where the second atom is located.
        res_id_2 (tuple): Tuple of (residue sequence number, insertion code) for the second atom.
        atom_name_2 (str): Name of the second atom.
        add_hydrogen_to_1 (bool): Add a hydrogen to atom_1 for distance calculation.
        add_hydrogen_to_2 (bool): Add a hydrogen to atom_2 for distance calculation.

    Returns:
        float: Distance between the specified atoms (or hydrogen atoms) in angstroms.
    """

    file_format = os.path.splitext(structure_file)[1][1:]

    if file_format.lower() == 'pdb':
        parser = PDBParser(QUIET=True)
    elif file_format.lower() == 'cif':
        parser = MMCIFParser(QUIET=True)
    else:
        raise ValueError("Unsupported file format. Use 'pdb' or 'cif'.")

    structure = parser.get_structure("protein", structure_file)

    # Locate the first atom
    chain_1 = structure[0][chain_id_1]
    chain_2 = structure[0][chain_id_2]

    residue_1 = find_residue_by_id(chain_1, res_id_1)
    residue_2 = find_residue_by_id(chain_2, res_id_2)

    atom_1 = get_atom_with_variations(residue_1, atom_name_1)
    atom_2 = get_atom_with_variations(residue_2, atom_name_2)
    print(f"Measuring bond distance between {chain_1}{residue_1}{atom_1} and {chain_2}{residue_2}{atom_2}.")

    # Determine coordinates for atom_1 and atom_2
    if add_hydrogen_to_1:
        neighbors_1 = get_covalent_neighbors(atom_1, residue_1)
        coord_1 = calculate_hydrogen_position(atom_1, neighbors_1)
    else:
        coord_1 = atom_1.coord

    if add_hydrogen_to_2:
        neighbors_2 = get_covalent_neighbors(atom_2, residue_2)
        coord_2 = calculate_hydrogen_position(atom_2, neighbors_2)
    else:
        coord_2 = atom_2.coord

    # Calculate the distance
    distance = np.linalg.norm(coord_1 - coord_2)
    return distance


class BondData(ZSData):
    def __init__(
        self,
        input_csv: str,
        scale_fit: str = "not_scaled",
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        mut_col_name: str = "mut",
        pos_col_name: str = "pos",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        seq_dir: str = "data/seq",
        zs_dir: str = "zs",
        struct_dir: str = "af3_joint/struct/PfTrpB-4bromo",
        # triad_dir: str = "triad",
        # triad_mut_dir: str = "mut_file",
        # triad_struct_dir: str = "struct_file",
        # withsub: bool = True,
        # chain_id: str = "A",
    ):

        super().__init__(
            input_csv=input_csv,
            scale_fit=scale_fit,
            combo_col_name=combo_col_name,
            var_col_name=var_col_name,
            mut_col_name=mut_col_name,
            pos_col_name=pos_col_name,
            seq_col_name=seq_col_name,
            fit_col_name=fit_col_name,
            # withsub=withsub,
            seq_dir=seq_dir,
            zs_dir=zs_dir,
        )

   