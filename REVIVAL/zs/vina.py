"""
A script for get docking scores from chai structures

must use vina conda env
"""

from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Union

import logging
import subprocess

import warnings

import os
import re

import shutil
from glob import glob
from tqdm import tqdm
from copy import deepcopy

import numpy as np
import pandas as pd

from rdkit.Chem.MolStandardize.rdMolStandardize import Uncharger
from pdbfixer import PDBFixer
from openmm.app import PDBFile, PDBxFile
from Bio import PDB
from Bio.PDB import PDBParser, PDBIO, Select, MMCIFIO
from Bio.PDB.Atom import Atom

from rdkit import Chem
from rdkit.Chem import AllChem

from MDAnalysis import Universe

from REVIVAL.global_param import LIB_INFO_DICT, ENZYME_INFO_DICT, AA_DICT
from REVIVAL.preprocess import ZSData
from REVIVAL.chem_helper import protonate_smiles, apply_mutation
from REVIVAL.util import (
    checkNgen_folder,
    get_file_name,
    get_protein_structure,
    get_chain_ids,
    get_chain_structure,
    calculate_chain_centroid,
    calculate_ligand_centroid,
    replace_residue_names_auto,
    save_hetatm_only
)


warnings.filterwarnings("ignore")


LIGAND_IONS = ["NA", "FE"]


def format_pdb(file_path: str, protein_dir: str, pH: float, regen=False):

    """
    Check if we have a structure,
    Clean this guy for docking and other random shit.
    """

    # First check if they passed a filename or the sequence name
    assert os.path.isfile(file_path), f"f{file_path} not a file"

    name = get_file_name(file_path)

    this_protein_dir = checkNgen_folder(os.path.join(protein_dir, name))

    protein_pdb_file = os.path.join(this_protein_dir, name + ".pdb")
    protein_pdbqt_file = os.path.join(this_protein_dir, name + ".pdbqt")

    if not os.path.isfile(protein_pdbqt_file) or regen:

        # remove and add back Hs for the protein only
        clean_one_pdb(file_path, protein_pdb_file, keep_chain="A")

        # Convert to pdbqt for vina
        pdb_to_pdbqt_protein(
            input_path=protein_pdb_file, output_path=protein_pdbqt_file, pH=pH
        )

    return protein_pdb_file, protein_pdbqt_file


##### helper functions for cleaning pdb #####


def save_clean_protein(s, toFile: str, keep_chain="A", keep_all_protein_chains=True):
    """
    First remove everything before adding everything back in.
    """

    class MySelect(Select):
        def accept_residue(self, residue, keep_chain=keep_chain):
            pdb, _, chain, (hetero, resid, insertion) = residue.full_id
            if keep_all_protein_chains or (chain == keep_chain):
                # if hetero == ' ' or hetero == 'H_MSE':
                if hetero == "H_MSE":
                    residue.id = (" ", resid, insertion)
                    return True
                elif hetero == " ":
                    return True
                else:
                    return False
            else:
                return False

        def accept_atom(self, atom):
            # remove altloc atoms.
            return not atom.is_disordered() or atom.get_altloc() == "A"

    if toFile[-4:] == ".cif":
        io = MMCIFIO()
    elif toFile[-4:] == ".pdb":
        io = PDBIO()
    else:
        print("toFile should end with .cif or .pdb")
    io.set_structure(s)
    io.save(toFile, MySelect())


def remove_hydrogen_pdb(pdbFile: str, toFile: str) -> None:

    """
    Remove hydrogens from a pdb file.
    """

    s = get_protein_structure(pdbFile)

    class NoHydrogen(Select):
        def accept_atom(self, atom):
            if atom.element == "H" or atom.element == "D":
                return False
            return True

    io = MMCIFIO() if toFile[-4:] == ".cif" else PDBIO()
    io.set_structure(s)
    io.save(toFile, select=NoHydrogen())


def clean_one_pdb(proteinFile: str, toFile: str, keep_chain="keep_all"):

    """
    Clean and then remove stuff from a pdb file and then read-add in missing things.
    """

    s = get_protein_structure(proteinFile)

    if keep_chain == "keep_all":
        save_clean_protein(s, toFile, keep_all_protein_chains=True)
    else:
        save_clean_protein(
            s, toFile, keep_chain=keep_chain, keep_all_protein_chains=False
        )

    pdbFile = toFile
    fixed_pdbFile = toFile
    # Remove then readd hydrogens
    remove_hydrogen_pdb(pdbFile, fixed_pdbFile)
    fixer = PDBFixer(filename=fixed_pdbFile)
    fixer.removeHeterogens()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms(seed=0)

    try:
        fixer.addMissingHydrogens()
    except Exception as e:
        print(f"Skipping addMissingHydrogens due to error: {e}")

    if pdbFile[-3:] == "pdb":
        PDBFile.writeFile(
            fixer.topology, fixer.positions, open(fixed_pdbFile, "w"), keepIds=True
        )
    elif pdbFile[-3:] == "cif":
        PDBxFile.writeFile(
            fixer.topology, fixer.positions, open(fixed_pdbFile, "w"), keepIds=True
        )
    else:
        raise "protein is not pdb or cif"


def pdb_to_pdbqt_protein(input_path: str, output_path=None, pH: float = 7.4):

    """
    Convert a pdb file to a pdbqt file.
    """

    # Need to first remove stuff that is sometimes added by
    lines = []
    with open(input_path, "r+") as fin:
        for line in fin:
            if (
                line.split(" ")[0] not in ["ENDBRANCH", "BRANCH", "ROOT", "ENDROOT"]
                and "Fe" not in line
            ):  # Add in the removal of the Iron bit
                lines.append(line)
    with open(input_path, "w+") as fout:
        for line in lines:
            fout.write(line)

    output_path = output_path if output_path else input_path.replace(".pdb", ".pdbqt")
    os.system(
        f"obabel {input_path} -xr -p {pH} --partialcharge gasteiger -O {output_path}"
    )
    # Now we also want to be cheeky and remove any secondary model parts from the file
    # This is a hacky way to keep a bound heme or something, seems to work fine.
    lines = []
    with open(output_path, "r+") as fin:
        for line in fin:
            if line.split(" ")[0] not in ["MODEL", "TER", "ENDMDL", "REMARK"]:
                lines.append(line)
    with open(output_path, "w+") as fout:
        for line in lines:
            if "ENDMDL" not in line:
                fout.write(line)
        fout.write("TER\n")


### format ligand ###


def prepare_ion(
    smiles: str, element: str, output_dir: str, input_struct_path: str, pH: float
) -> str:
    """Handle ion extraction and conversion to PDBQT."""

    if element is None or element == "":
        element = re.search(r"\[([A-Za-z]+)", smiles).group(1)

    ion_pdb_file = os.path.join(output_dir, f"{element}.pdb")
    ion_pdbqt_file = os.path.join(output_dir, f"{element}.pdbqt")

    # Extract ion from PDB
    extract_ions(
        input_file_path=input_struct_path, output_file_path=ion_pdb_file, ion=element
    )

    # Convert to PDBQT using Open Babel
    subprocess.run(
        ["obabel", ion_pdb_file, "-O", ion_pdbqt_file, "--metal"], check=True
    )

    return ion_pdbqt_file


def ligand_smiles2pdbqt(
    smiles: str, ligand_sdf_file: str, ligand_pdbqt_file: str, pH: float, regen: bool = False
) -> str:
    """Generate 3D coordinates, protonate, and save ligand."""
    if not regen and os.path.isfile(ligand_pdbqt_file):
        print(f"Skipping ligand {ligand_pdbqt_file} as it already exists and regen false.")
        return

    smiles = Chem.CanonSmiles(smiles, useChiral=True)
    protonated_smiles = protonate_smiles(smiles, pH=pH)
    mol = Chem.MolFromSmiles(protonated_smiles)
    # uncharger = Uncharger()
    # mol = uncharger.uncharge(mol)
    mol = Chem.AddHs(mol)

    AllChem.EmbedMolecule(mol)
    writer = Chem.SDWriter(ligand_sdf_file)
    writer.write(mol)
    writer.close()

    subprocess.run(["obabel", ligand_sdf_file, "-O", ligand_pdbqt_file], check=True)

    return ligand_pdbqt_file


def format_ligand(
    smiles: str,
    ligand_name: str,
    var_dir: str,
    input_struct_path: str,
    pH: float = 7.4,
    from_pdb: bool = True,
    if_substrate: bool = True,
    ligand_chain_id: str = None,
    ligand_info: list = None,
    substruct_criteria: str = None,
    regen: bool = False,
    target_list_addh: list = None,
    double_bond_pairs=None,
    branch_info=None,
) -> str:

    """
    Format a ligand for docking, supporting both ions and standard molecules.
    TODO need to extract the pdbqt file for the ligand directly from the pdb file
    unless otherwise specified.

    take in the chain(s) for the ligands,
    if was docked jointly will need to tease apart different ligands first

    if smiles_opt is not None, then the substrate will be converted from SMILES string

    Args:
        smiles (str): SMILES string of the ligand.
        ligand_name (str): Name of the ligand.
        ligand_dir (str): Directory to save the prepared files.
        pH (float): pH for protonation.
        from_pdb (bool): Extract the ligand from PDB if True.
        input_struct_path (str): Path to the PDB file.
        regen (bool): Whether to regenerate files if they already exist.

    Returns:
        str: Paths to the PDBQT file
    """

    this_ligand_dir = checkNgen_folder(os.path.join(var_dir, ligand_name))

    ligand_pdbqt_file = os.path.join(this_ligand_dir, f"{ligand_name}.pdbqt")
    ligand_pdb_file = os.path.join(this_ligand_dir, f"{ligand_name}.pdb")
    ligand_pdb_temp_file = os.path.join(this_ligand_dir, f"{ligand_name}_temp.pdb")
    ligand_pdb_prehydrogen_file = os.path.join(
        this_ligand_dir, f"{ligand_name}_prehydrogen.pdb"
    )
    ligand_sdf_file = os.path.join(this_ligand_dir, f"{ligand_name}.sdf")

    # Skip processing if file exists and regen is False
    if os.path.isfile(ligand_pdbqt_file) and not regen:
        print(
            f"Skipping ligand {ligand_name} as {ligand_pdbqt_file} already exists and regen false."
        )
        return ligand_pdbqt_file

    print(f"Processing ligand {ligand_name}")

    # Special handling for ions
    print("Checking for ions...")
    simple_ion = ligand_name.upper()[:2]
    if simple_ion in LIGAND_IONS or (smiles and smiles[0] == "[" and smiles[-1] == "]"):
        print(f"Processing ion {ligand_name}")

        return prepare_ion(
            smiles=smiles,
            element=simple_ion if simple_ion in LIGAND_IONS else None,
            output_dir=this_ligand_dir,
            input_struct_path=input_struct_path,
            pH=pH,
        )

    # try to directly convert from pdb
    if from_pdb:
        print(f"Processing ligand {ligand_name} from {input_struct_path}")

        if ligand_info is not None:

            # Extract ligand from PDB file
            extract_substruct(
                input_file_path=os.path.abspath(input_struct_path),
                output_substruct_path=os.path.abspath(ligand_pdb_temp_file),
                info_list=ligand_info,
                chain_id=ligand_chain_id,
                criteria=substruct_criteria,  # only include the ligand
            )

        else:
            get_chain_structure(
                input_file_path=input_struct_path,
                output_file_path=ligand_pdb_temp_file,
                chain_id=ligand_chain_id,
            )

        # get rid of ions as they will be handled separately
        remove_ions(
            input_file_path=ligand_pdb_temp_file,
            output_file_path=ligand_pdb_prehydrogen_file,
        )

        if target_list_addh is not None and "borane" not in ligand_name.lower():
            if "silane" in ligand_name.lower():
                cutoff = 2
            else:
                cutoff = 1.6
            add_hydrogens_to_atoms(
                input_pdb=ligand_pdb_prehydrogen_file,
                output_pdb=ligand_pdb_file,
                target_list=target_list_addh,
                bond_length=1.0,
                neighbor_cutoff=cutoff,
                double_bond_pairs=double_bond_pairs,
            )
            os.remove(ligand_pdb_prehydrogen_file)
        elif "borane" in ligand_name.lower():
            add_hydrogens_to_boron(
                input_pdb=ligand_pdb_prehydrogen_file,
                output_pdb=ligand_pdb_file,
                target_atom=target_list_addh[0],
                bond_length=1.0,
                neighbor_cutoff=1.6,
                double_bond_pairs=double_bond_pairs,
            )
        else:
            # rename the file
            os.rename(ligand_pdb_prehydrogen_file, ligand_pdb_file)

        # now remove the temp files
        os.remove(ligand_pdb_temp_file)

        # Construct command
        try:

            ligand_pdbqt_temp_file = os.path.join(
                this_ligand_dir, f"{ligand_name}_temp.pdbqt"
            )
            # make -xr optional depends on the ligand

            command = [
                "obabel",
                os.path.abspath(ligand_pdb_file),
                # "--assignatoms",
                # "--addh",
                # "--addexplicit",
                "-xr",
                # "-p", str(pH),
                "--partialcharge",
                "gasteiger",
                "-O",
                os.path.abspath(ligand_pdbqt_temp_file),
            ]

            if if_substrate and ("borane" not in ligand_name.lower()):
                command.remove("-xr")

            if ("borane" in ligand_name.lower()) or ("silane" in ligand_name.lower()):
                command.remove("--partialcharge")
                command.remove("gasteiger")

            # Execute command and capture output
            subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True,  # Raises CalledProcessError if the command fails
            )

            # now clean up
            if "borane" in ligand_name.lower():
                # remove the -H from the boron
                clean_boron_pdbqt_file(
                    ligand_pdbqt_temp_file,
                    ligand_pdbqt_file,
                    branch_B=branch_info["B"],
                    branch_C=branch_info["C"],
                )
            elif "silane" in ligand_name.lower():
                os.rename(ligand_pdbqt_temp_file, ligand_pdbqt_file)
            else:
                clean_pdbqt_file(ligand_pdbqt_temp_file, ligand_pdbqt_file)

        except subprocess.CalledProcessError as e:
            # Handle errors
            print(f"Error during conversion:\n{e.stderr}")
            return False

        except Exception as e:
            # Handle other exceptions
            print(f"An unexpected error occurred: {e}")
            return False

        # Validate output
        if (
            not os.path.isfile(ligand_pdbqt_file)
            or os.path.getsize(ligand_pdbqt_file) == 0
        ):
            raise ValueError(f"PDBQT file for {ligand_name} is empty or missing.")

        return ligand_pdbqt_file

    # if not from pdb, then we need to convert from smiles
    try:
        # Standard ligand processing
        ligand_smiles2pdbqt(
            smiles=smiles,
            ligand_sdf_file=ligand_sdf_file,
            ligand_pdbqt_file=ligand_pdbqt_file,
            pH=pH,
        )

    except Exception as e:
        raise ValueError(f"Error processing ligand '{ligand_name}': {e}")

    # Validate output
    if not os.path.isfile(ligand_pdbqt_file) or os.path.getsize(ligand_pdbqt_file) == 0:
        raise ValueError(f"PDBQT file for {ligand_name} is empty or missing.")

    return ligand_pdbqt_file


### helper functions for prep ligand ###


def rotation_matrix(axis, theta):
    """
    Return the 3x3 rotation matrix for rotating 'theta' radians about 'axis'.
    Uses Rodrigues' rotation formula.
    """
    axis = np.array(axis, dtype=float)
    axis /= np.linalg.norm(axis)
    x, y, z = axis
    c = np.cos(theta)
    s = np.sin(theta)
    C = 1 - c
    return np.array(
        [
            [x * x * C + c, x * y * C - z * s, x * z * C + y * s],
            [y * x * C + z * s, y * y * C + c, y * z * C - x * s],
            [z * x * C - y * s, z * y * C + x * s, z * z * C + c],
        ]
    )


def calculate_hydrogen_position(
    atom, neighbors, bond_length=1.0, double_bond_pairs=None
):
    """
    Calculate a new hydrogen coordinate for 'atom' based on:
      - number of existing neighbors
      - known double bonds (double_bond_pairs)
      - special case: O with 1 neighbor => offset angle for a more realistic OH geometry

    Args:
        atom (Bio.PDB.Atom.Atom): The central atom (will receive H).
        neighbors (list[Bio.PDB.Atom.Atom]): already-bonded neighbors of 'atom'.
        bond_length (float): approx distance for the new H–X bond.
        double_bond_pairs (dict or set or None):
            e.g. {("C3","N2"), ("N2","C3")}, marking certain bonds as double.

    Returns:
        np.ndarray: [x, y, z] of the newly placed hydrogen atom.
    """
    if double_bond_pairs is None:
        double_bond_pairs = {}

    element = atom.element.strip().upper()
    atom_name = atom.get_name().strip()
    atom_coord = np.array(atom.coord)
    n_neighbors = len(neighbors)

    # ----- Special case #1: O with 1 neighbor => "OH" with a small offset angle
    if element == "O" and n_neighbors == 1:
        neighbor_coord = np.array(neighbors[0].coord)

        # Vector from O to the neighbor
        vec = atom_coord - neighbor_coord
        norm_v = np.linalg.norm(vec)
        if norm_v < 1e-12:
            # degenerate fallback
            vec = np.array([1.0, 0.0, 0.0])
        else:
            vec /= norm_v

        # We'll rotate 'vec' by ~70° so the final angle is ~110° from the neighbor
        # This is a simple approximation for an -OH bond angle.
        angle_deg = 70.0
        angle_rad = np.radians(angle_deg)

        # Find a perpendicular axis to rotate around
        axis = np.cross(vec, [0, 0, 1])
        if np.linalg.norm(axis) < 1e-12:
            axis = np.cross(vec, [0, 1, 0])

        # Create rotation matrix
        R = rotation_matrix(axis, angle_rad)
        # Rotate the vector
        rotated_vec = R.dot(vec)
        # Scale by bond length
        H_coord = atom_coord + rotated_vec * bond_length
        return H_coord

    # ----- Special case #2: C with 2 neighbors + double bond => sp2 forced
    if element == "C" and n_neighbors == 2:
        neighbor_names = [n.get_name().strip() for n in neighbors]
        is_double_bonded = False
        for nbr_name in neighbor_names:
            if (atom_name, nbr_name) in double_bond_pairs:
                is_double_bonded = True
                break

        if is_double_bonded:
            # sp2 approach: ~120° from the two neighbors
            v1 = np.array(neighbors[0].coord) - atom_coord
            v2 = np.array(neighbors[1].coord) - atom_coord
            v1 /= np.linalg.norm(v1)
            v2 /= np.linalg.norm(v2)
            direction = -(v1 + v2)
            norm_d = np.linalg.norm(direction)
            if norm_d < 1e-12:
                # fallback
                direction = np.cross(v1, [0, 0, 1])
                if np.linalg.norm(direction) < 1e-12:
                    direction = np.array([1, 0, 0])
            else:
                direction /= norm_d
            return atom_coord + direction * bond_length

    # ----- Fallback: sp / sp2 / sp3 standard
    if n_neighbors == 1:
        # sp (linear)
        vec = atom_coord - np.array(neighbors[0].coord)
        nrm = np.linalg.norm(vec)
        if nrm < 1e-12:
            vec = np.array([1, 0, 0])
        else:
            vec /= nrm
        return atom_coord + vec * bond_length

    elif n_neighbors == 2:
        # sp2 (planar)
        coords = [np.array(n.coord) for n in neighbors]
        v1 = coords[0] - atom_coord
        v2 = coords[1] - atom_coord
        v1 /= np.linalg.norm(v1)
        v2 /= np.linalg.norm(v2)
        direction = -(v1 + v2)
        norm_d = np.linalg.norm(direction)
        if norm_d < 1e-12:
            direction = np.cross(v1, [0, 0, 1])
            if np.linalg.norm(direction) < 1e-12:
                direction = np.array([1, 0, 0])
        else:
            direction /= norm_d
        return atom_coord + direction * bond_length

    elif n_neighbors == 3:
        # sp3 (tetrahedral)
        coords = [np.array(n.coord) for n in neighbors]
        centroid = np.mean(coords, axis=0)
        direction = atom_coord - centroid
        norm_d = np.linalg.norm(direction)
        if norm_d < 1e-12:
            direction = np.array([1, 0, 0])
        else:
            direction /= norm_d
        return atom_coord + direction * bond_length

    # else no geometry known
    raise ValueError(
        f"Unsupported geometry for atom '{atom_name}' with {n_neighbors} neighbors."
    )


def calculate_hydrogen_positions_sp3(atom, neighbors, bond_length=1.0):
    """
    Calculate hydrogen positions for sp3 hybridization with tetrahedral geometry.

    Args:
        atom: The central atom (e.g., B in your case).
        neighbors: List of neighboring atoms already bonded to the central atom.
        bond_length: Desired bond length for the new H–X bond.

    Returns:
        List of np.ndarray: Coordinates for the new hydrogens.
    """
    atom_coord = np.array(atom.coord)
    neighbor_coords = [np.array(neighbor.coord) for neighbor in neighbors]

    # Case: 1 neighbor (sp³ hybridized, tetrahedral geometry for 3 Hs)
    if len(neighbor_coords) == 1:
        neighbor_vec = neighbor_coords[0] - atom_coord
        neighbor_vec /= np.linalg.norm(neighbor_vec)

        # Generate two orthogonal vectors to the neighbor_vec
        ortho_vec1 = np.cross(neighbor_vec, [1, 0, 0])
        if (
            np.linalg.norm(ortho_vec1) < 1e-6
        ):  # Handle edge case where neighbor_vec || [1, 0, 0]
            ortho_vec1 = np.cross(neighbor_vec, [0, 1, 0])
        ortho_vec1 /= np.linalg.norm(ortho_vec1)

        ortho_vec2 = np.cross(neighbor_vec, ortho_vec1)
        ortho_vec2 /= np.linalg.norm(ortho_vec2)

        # Scale the bond length to place hydrogens correctly
        tetrahedral_dirs = [
            (-neighbor_vec + ortho_vec1 + ortho_vec2),  # First hydrogen direction
            (-neighbor_vec - ortho_vec1 + ortho_vec2),  # Second hydrogen direction
            (-neighbor_vec - ortho_vec2),  # Third hydrogen direction
        ]

        # Normalize and scale directions to bond length
        h_positions = [
            atom_coord + bond_length * (direction / np.linalg.norm(direction))
            for direction in tetrahedral_dirs
        ]
        return h_positions

    # Case: Unsupported geometries
    else:
        raise ValueError(
            f"Unsupported geometry for atom with {len(neighbor_coords)} neighbors."
        )


def add_hydrogens_to_boron(
    input_pdb: str,
    output_pdb: str,
    target_atom,
    double_bond_pairs,
    bond_length=1.1,
    neighbor_cutoff=1.6,
):
    """
    Adds three hydrogens to a specified boron atom in the input PDB file.
    Args:
        input_pdb: Path to the input PDB file.
        output_pdb: Path to save the modified PDB file.
        target_atom: Tuple (chain_id, residue_name, atom_name) identifying the target boron atom.
        bond_length: Distance for the H-B bond.
        neighbor_cutoff: Cutoff distance to identify neighbors.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("my_structure", input_pdb)
    model = structure[0]

    chain_id, residue_name, atom_name = target_atom

    # Locate the target atom
    chain = model[chain_id]
    residue = next(res for res in chain if res.resname.strip() == residue_name)
    atom = next(a for a in residue if a.get_name().strip() == atom_name)

    # Find neighbors of the target atom
    neighbors = []
    for neighbor in residue.get_atoms():
        if neighbor is atom:
            continue
        dist = atom - neighbor
        if dist < neighbor_cutoff:
            neighbors.append(neighbor)

    # Calculate new hydrogen positions
    hydrogen_positions = calculate_hydrogen_positions_sp3(atom, neighbors, bond_length)

    # Add hydrogens to the structure
    for idx, h_coord in enumerate(hydrogen_positions, start=1):
        h_name = f"H{idx}"[:4].ljust(4)
        new_atom = Atom(
            h_name.strip(),
            h_coord,
            bfactor=0.0,
            occupancy=1.0,
            altloc=" ",
            fullname=h_name,
            serial_number=len(residue) + idx,
            element="H",
        )
        residue.add(new_atom)
        print(f"Added hydrogen {h_name.strip()} to {atom_name}")

    # Save the modified structure
    io = PDBIO()
    io.set_structure(structure)

    write_pdb_with_doublebonds(structure, output_pdb, double_bond_pairs)

    # io.save(output_pdb)
    print(f"Saved modified PDB to: {output_pdb}")


def add_hydrogens_to_atoms(
    input_pdb: str,
    output_pdb: str,
    target_list,
    bond_length: float = 1.0,
    neighbor_cutoff: float = 1.6,
    double_bond_pairs=None,
):
    """
    1) Reads the PDB
    2) Finds the largest existing H index
    3) For each (chain_id, residue_name, atom_name) in target_list:
         - Locates that atom
         - Finds neighbors
         - Calls 'calculate_hydrogen_position' (which uses double_bond_pairs if needed)
         - Adds one new hydrogen
    4) Writes out a custom PDB that includes only the CONECT lines for the pairs in double_bond_pairs.
       (No distance-based guess for other bonds).

    We assume double_bond_pairs is a list of 2-tuples:
      [ (("B","LIG","C3"), ("B","LIG","N2")), ... ]
    Each element is a pair of (chain, resname, atom_name).
    """
    if double_bond_pairs is None:
        double_bond_pairs = {}

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("my_structure", input_pdb)
    model = structure[0]

    # --- Step A: find the largest existing H index ---
    largest_h_num = 0
    pattern = re.compile(r"^H(\d+)$")
    for chain in model:
        for residue in chain:
            for atom in residue:
                if atom.element.strip().upper() == "H":
                    match_obj = pattern.match(atom.get_name().strip())
                    if match_obj:
                        val = int(match_obj.group(1))
                        if val > largest_h_num:
                            largest_h_num = val
    new_h_serial = largest_h_num + 1

    # --- Step B: For each target, add an H
    for (chain_id, residue_name, atom_name) in target_list:
        # locate chain
        if chain_id not in model:
            raise ValueError(f"Chain '{chain_id}' not found in PDB.")

        chain = model[chain_id]

        # find residue
        tgt_res = None
        for res in chain:
            if res.resname.strip() == residue_name:
                tgt_res = res
                break
        if not tgt_res:
            raise ValueError(
                f"Residue '{residue_name}' not found in chain '{chain_id}'."
            )

        # find atom
        tgt_atom = None
        for a in tgt_res:
            if a.get_name().strip() == atom_name:
                tgt_atom = a
                break
        if not tgt_atom:
            raise ValueError(
                f"Atom '{atom_name}' not in residue '{residue_name}' (chain '{chain_id}')"
            )

        # neighbors
        neighbors = []
        for a in tgt_res.get_atoms():
            if a is tgt_atom:
                continue
            dist = tgt_atom - a
            if dist < neighbor_cutoff:
                neighbors.append(a)

        # compute new H coordinate (using double_bond_pairs in the geometry if needed)
        new_coord = calculate_hydrogen_position(
            atom=tgt_atom,
            neighbors=neighbors,
            bond_length=bond_length,
            double_bond_pairs=double_bond_pairs,
        )

        # create a unique H name
        new_h_name = f"H{new_h_serial}"[:4].ljust(4)
        new_atom = Atom(
            new_h_name.strip(),
            new_coord,
            bfactor=0.0,
            occupancy=1.0,
            altloc=" ",
            fullname=new_h_name,
            serial_number=new_h_serial,
            element="H",
        )
        new_h_serial += 1

        tgt_res.add(new_atom)
        print(
            f"Added H '{new_h_name.strip()}' to {atom_name} in {residue_name} (chain {chain_id})"
        )

    if double_bond_pairs is not None:
        # --- Step C: Write updated PDB with specified double-bond CONECT lines
        write_pdb_with_doublebonds(
            structure=structure,
            output_pdb=output_pdb,
            double_bond_pairs=double_bond_pairs,
        )
    # if "silane" in input_pdb:
    #     temp_pdb = output_pdb.replace(".pdb", "_presifix.pdb")

    #     io = PDBIO()
    #     io.set_structure(structure)
    #     io.save(temp_pdb)

    #     update_si_type(input_pdb=temp_pdb, output_pdb=output_pdb)
    #     # os.remove(temp_pdb)

    else:
        # --- Step C: Write updated PDB
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb)

    print(f"\nSaved updated PDB to: {output_pdb}")


def write_pdb_with_doublebonds(structure, output_pdb, double_bond_pairs):
    """
    Write a PDB that includes ATOM/HETATM lines for all atoms in 'structure',
    plus CONECT lines ONLY for the pairs in 'double_bond_pairs'.

    Format for double_bond_pairs is a list of 2-tuples:
      [ (("B","LIG","C3"), ("B","LIG","N2")), (("B","LIG","C8"), ("B","LIG","N5")) ... ]

    We'll map each atom -> a unique serial, then write:
       HETATM (or ATOM) lines
       CONECT lines only for pairs we find in double_bond_pairs
    """
    # 1) Flatten all atoms. We'll store a key => (chain, resname, atomname)
    #    plus a big dictionary to get the assigned "serial"
    all_atoms = []
    atom_map = {}  # maps (chain, resname, atom_name) -> (atom_obj, serial)
    serial_counter = 1

    for chain in structure.get_chains():
        chain_id = chain.id
        for residue in chain:
            resname = residue.resname.strip()
            for atom in residue:
                a_name = atom.get_name().strip()
                key = (chain_id, resname, a_name)
                atom_map[key] = (atom, serial_counter)
                all_atoms.append(key)
                serial_counter += 1

    # 2) Write out the lines
    with open(output_pdb, "w") as outf:
        # Write HETATM lines
        for key in all_atoms:
            (atom_obj, serial) = atom_map[key]
            chain_id, resname, a_name = key
            x, y, z = atom_obj.coord
            element = atom_obj.element.strip().upper() or a_name[0]

            record_type = "HETATM"  # or "ATOM " if you prefer

            # We'll assume residue id is an integer. For real code, check insertion codes etc.
            # Biopython residue .id -> (hetero_flag, seq_number, insertion_code)
            # We'll guess no insertion code for simplicity
            seqnum = residue.id[1] if hasattr(atom_obj.parent, "id") else 1

            line = (
                f"{record_type:6s}{serial:5d} {a_name:<4s}"
                f"{resname:>4s} {chain_id:1s}"
                f"{seqnum:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}"
                f"{1.00:6.2f}{0.00:6.2f}          "
                f"{element:>2s}\n"
            )
            outf.write(line)

        # 3) For each double-bond pair, look up serials and write CONECT
        for pair in double_bond_pairs:
            # pair is e.g. (("B","LIG","C3"), ("B","LIG","N2"))
            a1, a2 = pair  # each is (chain_id, resname, atom_name)
            if a1 in atom_map and a2 in atom_map:
                s1 = atom_map[a1][1]
                s2 = atom_map[a2][1]
                outf.write(f"CONECT{s1:5d}{s2:5d}\n")
                # optionally also the reverse if you want symmetrical lines:
                outf.write(f"CONECT{s2:5d}{s1:5d}\n")
            else:
                print(
                    f"Warning: Could not find one of {a1} or {a2} in structure for CONECT."
                )

        outf.write("END\n")


def clean_boron_pdbqt_file(
    input_pdbqt: str, output_pdbqt: str, branch_B="B1", branch_C="C1"
):

    """
    Group B1 and its hydrogens into a branch structure, and group the rest with C1.

    Args:
        input_pdbqt (str): Path to the input PDBQT file.
        output_pdbqt (str): Path to the output PDBQT file.
    """
    try:
        with open(input_pdbqt, "r") as infile, open(output_pdbqt, "w") as outfile:
            lines = infile.readlines()
            root_atoms = []
            branch_atoms = []
            branch_hydrogens = []

            b1_serial = None
            c1_serial = None

            # Process each line and categorize atoms
            for line in lines:
                if line.startswith("TER"):
                    continue  # Skip TER lines
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    atom_serial = int(line[6:11].strip())
                    atom_name = line[12:16].strip()

                    if atom_name == branch_B:
                        # update B to C
                        atom_type = line[77:79].strip()
                        assert atom_type == "B"  # Replace 'B' with 'C'
                        line = f"{line[:77]}C {line[79:]}"
                        branch_atoms.append(line)
                        b1_serial = atom_serial
                    elif atom_name == branch_C:
                        root_atoms.append(line)
                        c1_serial = atom_serial
                    elif line[77:79].strip() == "HD":
                        branch_hydrogens.append(line)
                    else:
                        root_atoms.append(line)

            if b1_serial is None or c1_serial is None:
                raise ValueError("Could not find B1 or C1 in the input PDBQT file.")

            # Write header
            for line in lines:
                if (
                    not line.startswith("ATOM")
                    and not line.startswith("HETATM")
                    and not line.startswith("TER")
                ):
                    outfile.write(line)

            # Write ROOT section
            outfile.write("ROOT\n")
            for atom in root_atoms:
                outfile.write(atom)
            outfile.write("ENDROOT\n")

            # Write BRANCH section
            outfile.write(f"BRANCH   {c1_serial}   {b1_serial}\n")
            for atom in branch_atoms:
                outfile.write(atom)
            for hydrogen in branch_hydrogens:
                outfile.write(hydrogen)
            outfile.write(f"ENDBRANCH   {c1_serial}   {b1_serial}\n")

            # Add torsional degrees of freedom
            outfile.write("TORSDOF 1\n")

        print(f"Processed PDBQT file saved to: {output_pdbqt}")

    except Exception as e:
        print(f"Error processing file: {e}")


def clean_pdbqt_file(input_pdbqt, output_pdbqt):
    """
    Cleans up a PDBQT file:
    - Replaces 'ACE A   0' with 'UNL     1' while maintaining proper formatting.
    - Ensures spacing is correct for all fields.
    """
    with open(input_pdbqt, "r") as infile, open(output_pdbqt, "w") as outfile:
        for line in infile:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Check and replace the residue name and residue number
                resname = line[17:20].strip()
                chain_id = line[21].strip()
                resseq = line[22:26].strip()

                if resname == "ACE" and chain_id == "A" and resseq == "0":
                    # Replace with 'UNL     1'
                    line = f"{line[:17]}UNL     1{line[26:]}"

            # Write the modified (or unmodified) line to the output
            outfile.write(line)

    print(f"Cleaned PDBQT file saved to: {output_pdbqt}")


def remove_ions(input_file_path: str, output_file_path: str) -> None:
    """
    Remove lines containing ions from a PDB or CIF file.

    Args:
        input_file_path (str): Path to the input PDB or CIF file.
        output_file_path (str): Path to save the file without ions.
    """
    print(f"Removing ions from {input_file_path}")

    with open(input_file_path, "r") as infile, open(output_file_path, "w") as outfile:
        for line in infile:
            if not (line.startswith("ATOM") or line.startswith("HETATM")) or not any(
                ion in line[76:78].strip() for ion in LIGAND_IONS
            ):
                outfile.write(line)


def extract_ions(input_file_path: str, output_file_path: str, ion: str) -> None:
    """
    Extract lines starting with ATOM containing the specified ion from a PDB or CIF file.

    Args:
        input_file_path (str): Path to the input PDB or CIF file.
        output_file_path (str): Path to save the extracted ion lines.
        ion (str): Name of the ion to extract (e.g., "Na", "Cl").
    """

    ion_found = False

    with open(input_file_path, "r") as infile, open(output_file_path, "w") as outfile:
        for line in infile:
            if line.startswith("ATOM") or line.startswith("HETATM"):

                fields = line.split()
                if len(fields) >= 4:
                    atom_name = fields[2].strip()  # Atom name
                    residue_name = fields[3].strip()  # Residue name
                    element_name = fields[-1].strip()
                    if (
                        atom_name.upper() == ion.upper()
                        or residue_name.upper() == ion.upper()
                        or element_name.upper() == ion.upper()
                    ):
                        outfile.write(line)
                        ion_found = True

    if not ion_found:
        print(f"No ions matching '{ion}' were found in the file.")
    else:
        print(
            f"Ions matching '{ion}' were successfully extracted to {output_file_path}."
        )


def save_full_hierarchy_atoms(atoms_with_context, output_file_path):
    """
    Save a PDB file containing atoms with full hierarchy (model, chain, residue, atom).

    Args:
        atoms_with_context (list): List of tuples (model, chain, residue, atom).
        output_file_path (str): Path to save the filtered PDB file.
    """

    # Create a custom selector
    class FullHierarchySelector(Select):
        def __init__(self, atoms_with_context):
            self.selected_atoms = set(atom for _, _, _, atom in atoms_with_context)

        def accept_atom(self, atom):
            # Include only selected atoms
            return atom in self.selected_atoms

    # Use the structure from the first atom's model as a base
    if atoms_with_context:

        structure = atoms_with_context[0][0].get_parent()

        # Write the filtered structure
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_file_path, select=FullHierarchySelector(atoms_with_context))


def extract_substruct(
    input_file_path: str,
    output_substruct_path: str,
    info_list: list,
    chain_id: str,
    criteria: str = "include",
):
    """
    Extract substructure from a PDB or CIF file based on chain, resn, and atom name.

    Args:
        input_file_path (str): Path to the input PDB or CIF file.
        output_substruct_path (str): Path to save the extracted substrate structure.
        info_list (list of tuples): List of (chain_id, resn, name) for selection.
        criteria (str): "include" to keep specified atoms, "exclude" to remove them.
    """

    structure = get_protein_structure(input_file_path)

    substrate_atoms = []
    rest_atoms = []

    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    for atom in residue:
                        match = any(
                            chain.id == chain_id
                            and residue.resname == resn
                            and atom.name == name
                            for chain_id, resn, name in info_list
                        )
                        if match:
                            substrate_atoms.append((model, chain, residue, atom))
                        else:
                            rest_atoms.append((model, chain, residue, atom))

    print(f"Substrate atoms: {substrate_atoms}")
    print(f"Rest atoms: {rest_atoms}")

    if criteria == "include":
        save_full_hierarchy_atoms(substrate_atoms, output_substruct_path)
    else:
        save_full_hierarchy_atoms(rest_atoms, output_substruct_path)

    print(f"Substructure saved to {output_substruct_path}")


### docking ###
def dock(
    input_struct_path: str,
    dock_opt: str,  #  ie "substrate",
    score_only: bool,  # = True,
    lib_name: str = None,
    cofactor_dets: str = "cofactor",
    vina_dir: str = "zs/vina",
    residues4centriod: list = None,
    from_pdb: bool = True,
    pH: float = 7.4,
    size_x=20.0,
    size_y=20.0,
    size_z=20.0,
    num_modes=9,
    exhaustiveness=32,
    regen=False,
    rerun=False,
    num_cpus: int = None,
    seed=42,
):

    """
    Dock a substrate with multiple cofactors to a protein using vina.

    Args:
    - input_struct_path (str): Path to the input PDB file.
        ie. zs/chai/struct_joint/PfTrpB-4bromo/I165A:I183A:Y301V/I165A:I183A:Y301V_0.cif
        ie. zs/af3/struct_joint/PfTrpB-4bromo/i165a_i183a_y301v/seed-1_sample-0/model.cif
        ie. zs/af3/struct_joint/PfTrpB-4bromo/i165a_i183a_y301v/i165a_i183a_y301v_model.cif

    """

    if "joint" in input_struct_path:
        input_struct_dock_opt = "joint"
    elif "seperate" in input_struct_path:
        input_struct_dock_opt = "seperate"
    else:
        raise ValueError("Neither 'joint' nor 'seperate' found in clean_input_file")

    struct_tpye = "chai" if "chai" in input_struct_path else "af3"

    if struct_tpye == "chai":
        var_name = get_file_name(input_struct_path)
        if lib_name is None:
            lib_name = os.path.basename(
                os.path.dirname(os.path.dirname(input_struct_path))
            )

    elif struct_tpye == "af3":
        if lib_name is None:
            lib_name = input_struct_path.split(f"{input_struct_dock_opt}/")[-1].split(
                "/"
            )[0]
            print(lib_name)
        var_base_name = (
            input_struct_path.split(f"{lib_name}/")[-1]
            .split("/")[0]
            .replace("_", ":")
            .upper()
        )
        rep_numb = (
            input_struct_path.split("seed-1_sample-")[-1].split("/")[0]
            if "seed-1_sample-" in input_struct_path
            else "agg"
        )
        var_name = f"{var_base_name}_{rep_numb}"
        print(var_name)

    struct_dets = input_struct_path.split("zs/")[-1].split(lib_name)[0]
    var_dir = checkNgen_folder(os.path.join(vina_dir, struct_dets, lib_name, var_name))

    lib_info = LIB_INFO_DICT[lib_name]
    substrate_name = lib_info["substrate"]
    output_opt = "score_only" if score_only else "docked"

    # Auxiliary files
    vina_logfile = os.path.join(
        var_dir, f"{var_name}-{substrate_name}-{dock_opt}-{output_opt}_log.txt"
    )
    docked_ligand_pdb = os.path.join(
        var_dir, f"{var_name}-{substrate_name}-{dock_opt}-{output_opt}.pdb"
    )

    if not rerun and os.path.isfile(vina_logfile):
        print(f"{vina_logfile} exists and rerun is False. Skipping")
        return vina_logfile

    # convert cif to pdb
    pre_clean_input_file = os.path.join(var_dir, var_name + "_raw.pdb")

    # obabel input.cif -O output.pdb
    cmd = f"obabel {input_struct_path} -O {pre_clean_input_file}"
    subprocess.run(cmd, shell=True)

    # clean input file
    clean_input_file = os.path.join(var_dir, var_name + ".pdb")

    replace_residue_names_auto(
        input_file=pre_clean_input_file,
        output_file=clean_input_file,
    )

    # Process the protein
    protein_pdb_file, protein_pdbqt = format_pdb(
        file_path=clean_input_file, protein_dir=var_dir, pH=pH, regen=regen
    )

    if score_only:
        from_pdb = True

    ligand_chain_ids = get_chain_ids(pdb_file_path=clean_input_file)[1:]

    ligand_info = lib_info[f"{substrate_name}-info"]

    if len(ligand_chain_ids) == 1:
        ligand_chain_id = ligand_chain_ids[0]
        cofactor_chain_id = ligand_chain_id

    elif len(ligand_chain_ids) == 2:
        ligand_chain_id = ligand_info[0][0]
        cofactor_chain_id = [i for i in ligand_chain_ids if i != ligand_chain_id][0]
    else:
        raise ValueError("Substrate chains not implemented beyond 2")

    substrate_hinfo = lib_info.get(f"substrate-addH_{struct_tpye}", None)
    substrate_db_pairs = lib_info.get(
        f"substrate-double_bond_pairs_{struct_tpye}", None
    )
    branch_info = lib_info.get(f"substrate_branches_{struct_tpye}", None)

    # Process the main ligand
    ligand_pdbqt = format_ligand(
        smiles=lib_info["substrate-smiles"],
        ligand_name=substrate_name,
        var_dir=var_dir,
        pH=pH,
        substruct_criteria="include",
        from_pdb=from_pdb,
        input_struct_path=clean_input_file,
        ligand_info=ligand_info,
        ligand_chain_id=ligand_chain_id,
        target_list_addh=substrate_hinfo,
        double_bond_pairs=substrate_db_pairs,
        regen=regen,
        branch_info=branch_info,
    )

    # Process all cofactors from the input pdb
    cofactor_pdbqts = []
    for (cofactor_name, cofactor_smiles) in zip(
        LIB_INFO_DICT[lib_name][cofactor_dets],
        LIB_INFO_DICT[lib_name][f"{cofactor_dets}-smiles"],
    ):

        print(
            f"Processing cofactor {cofactor_name} {cofactor_smiles} using {clean_input_file}"
        )

        if f"carbene-info_{input_struct_dock_opt}" in lib_info:
            # need to get carbene and then heme and Fe
            carbene_info = lib_info[f"carbene-info_{input_struct_dock_opt}"]
            # TODO to check H add and pdbqt gen
            cofactor_pdbqts.append(
                format_ligand(
                    smiles=cofactor_smiles,
                    ligand_name="carbene",
                    var_dir=var_dir,
                    pH=pH,
                    substruct_criteria="include",
                    from_pdb=from_pdb,
                    input_struct_path=clean_input_file,
                    ligand_info=carbene_info,
                    ligand_chain_id=cofactor_chain_id,
                    regen=regen,
                )
            )

            if input_struct_dock_opt == "joint":
                heme_exclude_info = carbene_info + ligand_info
            else:
                heme_exclude_info = carbene_info

            cofactor_pdbqts.append(
                format_ligand(
                    smiles=cofactor_smiles,
                    ligand_name="heme",
                    var_dir=var_dir,
                    pH=pH,
                    substruct_criteria="exclude",
                    from_pdb=from_pdb,
                    input_struct_path=clean_input_file,
                    if_substrate=False,
                    ligand_info=heme_exclude_info,
                    ligand_chain_id=cofactor_chain_id,
                    regen=regen,
                )
            )

            # handle Fe separately
            cofactor_pdbqts.append(
                format_ligand(
                    smiles="[Fe2+]",
                    ligand_name="Fe",
                    var_dir=var_dir,
                    pH=pH,
                    from_pdb=from_pdb,
                    input_struct_path=clean_input_file,
                    ligand_info=None,
                    ligand_chain_id=cofactor_chain_id,
                    regen=regen,
                )
            )

        # trpbs
        else:
            cofactor_pdbqts.append(
                format_ligand(
                    smiles=cofactor_smiles,
                    ligand_name=cofactor_name,
                    var_dir=var_dir,
                    pH=pH,
                    substruct_criteria="exclude",
                    from_pdb=from_pdb,
                    input_struct_path=clean_input_file,
                    if_substrate=False,
                    ligand_chain_id=cofactor_chain_id,
                    ligand_info=ligand_info,
                    target_list_addh=lib_info[
                        f"cofactor-addH_{struct_tpye}_{input_struct_dock_opt}"
                    ],
                    double_bond_pairs=lib_info[
                        f"cofactor-double_bond_pairs_{struct_tpye}_{input_struct_dock_opt}"
                    ],
                    regen=regen,
                )
            )

    print(f"cofactor_pdbqts: {cofactor_pdbqts}")

    conf_path = os.path.join(
        var_dir, f"{var_name}-{substrate_name}-{dock_opt}_conf.txt"
    )

    if not os.path.exists(conf_path):

        # Create config file
        make_config_for_vina(
            input_pdb_path=input_struct_path,
            protein_pdbqt=protein_pdbqt,
            ligand_pdbqt=ligand_pdbqt,
            cofactor_pdbqts=cofactor_pdbqts,
            conf_path=conf_path,
            residues=residues4centriod,
            substrate_chain_ids=ligand_chain_ids,
            dock_opt=dock_opt,
            size_x=size_x,
            size_y=size_y,
            size_z=size_z,
            num_modes=num_modes,
            exhaustiveness=exhaustiveness,
        )

    # Construct the base command
    cmd_list = ["vina", "--config", conf_path, "--out", docked_ligand_pdb]

    if score_only:
        cmd_list += ["--score_only"]
    else:
        cmd_list += ["--seed", str(seed)]
        if num_cpus is not None:
            cmd_list += ["--cpu", str(num_cpus)]

    try:
        # Run the Vina command
        cmd_return = subprocess.run(
            cmd_list,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=True,  # Ensures subprocess raises an exception on error
        )

        # Write the Vina output to the log file
        with open(vina_logfile, "w") as fout:
            fout.write(cmd_return.stdout.decode("utf-8"))

        print(
            f"Vina run {'(score-only)' if score_only else ''} completed successfully. Log saved to: {vina_logfile}"
        )
        if not score_only:
            print(f"Docked ligand saved to: {docked_ligand_pdb}")

    except subprocess.CalledProcessError as e:
        print(f"Error during Vina run: {e}")
        with open(vina_logfile, "w") as fout:
            fout.write(e.stdout.decode("utf-8"))
        raise

    if not score_only:
        # Convert the docked PDBQT file to PDB if not in score-only mode
        convert_pdbqt_to_pdb(
            pdbqt_file=docked_ligand_pdb,
            pdb_file=docked_ligand_pdb,
            disable_bonding=True,
        )

    return vina_logfile


def dock_lib_parallel(
    struct_dir: str,
    dock_opt: str,  #  ie "substrate",
    score_only: bool,  # = True,
    cofactor_dets: str = "cofactor",
    vina_dir: str = "zs/vina",
    residues4centriod: list = None,
    from_pdb: bool = True,
    pH: float = 7.4,
    size_x=20.0,
    size_y=20.0,
    size_z=20.0,
    num_modes=9,
    exhaustiveness=32,
    regen=False,
    rerun=False,
    seed=42,
    num_cpus=None,  # for each dock function
    max_workers=24,  # Number of parallelized variants to be docked
):
    """
    A function to dock all generated chai structures and get scores in parallel.

    Args:
    - struct_dir (str): Path to the directory containing the structures.
        ie zs/chai/struct_joint/PfTrpB-4bromo
    """

    lib_name = os.path.basename(struct_dir)

    print(f"Docking {struct_dir} with {dock_opt}")

    struct_dets = struct_dir.split("zs/")[-1].split(lib_name)[0]

    checkNgen_folder(os.path.join(vina_dir, struct_dets, lib_name))

    if "af3" in struct_dir:
        var_path_agg_path = sorted(glob(os.path.join(struct_dir, "*", "*.cif")))
        var_rep_path = sorted(glob(os.path.join(struct_dir, "*", "*", "*.cif")))
        var_paths = var_path_agg_path + var_rep_path
    else:
        var_paths = sorted(glob(os.path.join(struct_dir, "*", "*.cif")))

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(
                dock,
                input_struct_path=var_path,
                dock_opt=dock_opt,
                score_only=score_only,
                lib_name=lib_name,
                cofactor_dets=cofactor_dets,
                vina_dir=vina_dir,
                residues4centriod=residues4centriod,
                from_pdb=from_pdb,
                pH=pH,
                size_x=size_x,
                size_y=size_y,
                size_z=size_z,
                num_modes=num_modes,
                exhaustiveness=exhaustiveness,
                regen=regen,
                rerun=rerun,
                num_cpus=num_cpus,
                seed=seed,
            )
            for var_path in var_paths
        ]

        for future in tqdm(as_completed(futures), total=len(var_paths)):
            try:
                future.result()  # Raises an exception if the task failed
            except Exception as e:
                print(e)


def make_config_for_vina(
    input_pdb_path: str,
    protein_pdbqt: str,
    ligand_pdbqt: str,
    cofactor_pdbqts: list,
    conf_path: str,
    residues: list = None,
    substrate_chain_ids: Union[list, str] = "B",
    dock_opt: str = "substrate",
    size_x=20.0,
    size_y=20.0,
    size_z=20.0,
    num_modes=9,
    exhaustiveness=32,
):
    """
    Create the config file for Vina, including cofactors if specified.

    Args:
    - input_pdb_path (str): Path to the input PDB file.
    - protein_pdbqt (str): Path to the protein PDBQT file.
        ie zs/vina/chai_joint/PfTrpB-4bromo/I165A:I183A:Y301V_0/I165A:I183A:Y301V_0/I165A:I183A:Y301V_0.pdbqt
    - ligand_pdbqt (str): Path to the ligand PDBQT file.
        ie zs/vina/chai_joint/PfTrpB-4bromo/I165A:I183A:Y301V_0/4bromo/4bromo.pdbqt
    - conf_path (str): Path to the output config file.
    - residues (list): List of residue positions to dock.
    - substrate_chain_ids (list): List of chain IDs for the substrate.
    - cofactor_files (list): List of paths to cofactor PDBQT files.
    - freeze_opt (str): Option for freezing cofactors,
        ie. cofactor means combining their pdbqt files with the enzyme pebqt file
            heme means just combining heme and Fe
    - size_x (float): Size of the docking box in the x dimension.
    - size_y (float): Size of the docking box in the y dimension.
    - size_z (float): Size of the docking box in the z dimension.
    - num_modes (int): Number of docking modes to generate.
    - exhaustiveness (int): Exhaustiveness of the docking
    """

    if residues is None:
        coords = calculate_chain_centroid(
            input_file=input_pdb_path, chain_ids=substrate_chain_ids
        )
    else:
        coords = []
        for position in residues:
            _, coordinates = get_coordinates_without_chain_id(
                input_pdb_path, int(position)
            )
            if coordinates is not None:
                coords.append(coordinates)
        coords = calculate_centroid(coords)

    var_dir = os.path.dirname(os.path.dirname(protein_pdbqt))
    var_name = get_file_name(protein_pdbqt)

    # treat all cofactors as ligands
    if dock_opt == "all":
        receptor_pdbqt = protein_pdbqt
        cofactor2dock = deepcopy(cofactor_pdbqts)

    # freeze cofactors as part of receptors
    elif dock_opt == "substrate":
        receptor_name = var_name + "_cofactors"
        # make an new subdir
        merge_receptor_dir = checkNgen_folder(os.path.join(var_dir, receptor_name))

        # Combine the protein and cofactor files
        receptor_pdbqt = os.path.join(merge_receptor_dir, f"{receptor_name}.pdbqt")
        receptor_pdb = os.path.join(merge_receptor_dir, f"{receptor_name}.pdb")

        merge_pdbqt(
            input_files=[protein_pdbqt] + cofactor_pdbqts,
            output_file_path=receptor_pdbqt,
        )

        # Convert the combined PDBQT file to PDB for downstream analysis
        convert_pdbqt_to_pdb(
            pdbqt_file=receptor_pdbqt, pdb_file=receptor_pdb, disable_bonding=True
        )

        cofactor2dock = None

    # freeze heme but dock carbene and the substrate
    elif dock_opt == "substrate+carbene":

        receptor_name = var_name + "_heme"
        # make an new subdir
        merge_receptor_dir = checkNgen_folder(os.path.join(var_dir, receptor_name))

        # Combine the protein and cofactor files
        receptor_pdbqt = os.path.join(merge_receptor_dir, f"{receptor_name}.pdbqt")
        receptor_pdb = os.path.join(merge_receptor_dir, f"{receptor_name}.pdb")

        cofactor2freeze = []

        for c in cofactor_pdbqts:
            if "carbene" in c.lower():
                cofactor2dock = [c]
            else:
                cofactor2freeze.append(c)

        print(f"cofactor2freeze: {cofactor2freeze}")
        print(f"cofactor2dock: {cofactor2dock}")

        # Combine the protein and cofactor files
        merge_pdbqt(
            input_files=[protein_pdbqt] + cofactor2freeze,
            output_file_path=receptor_pdbqt,
        )

        # Convert the combined PDBQT file to PDB
        convert_pdbqt_to_pdb(
            pdbqt_file=receptor_pdbqt, pdb_file=receptor_pdb, disable_bonding=True
        )

    else:
        raise ValueError(f"{dock_opt} invalid dock option")

    with open(conf_path, "w") as fout:
        fout.write(f"receptor = {receptor_pdbqt}\n")
        fout.write(f"ligand = {ligand_pdbqt}\n")

        # Include cofactors
        if cofactor2dock is not None:
            for cofactor_file in cofactor2dock:
                fout.write(f"ligand = {cofactor_file}\n")

        fout.write(f"center_x = {coords[0]}\n")
        fout.write(f"center_y = {coords[1]}\n")
        fout.write(f"center_z = {coords[2]}\n")
        fout.write(f"size_x = {size_x}\n")
        fout.write(f"size_y = {size_y}\n")
        fout.write(f"size_z = {size_z}\n")
        fout.write(f"num_modes = {num_modes}\n")
        fout.write(f"exhaustiveness = {exhaustiveness}\n")


def write_vina_affinities(protein_name, ligand_name, log_file_path, output_csv_path):
    # Read the content of the log file
    with open(log_file_path, "r") as file:
        log_content = file.readlines()

    # Parse the relevant lines
    data_lines = []
    add_lines = False
    for line in log_content:
        if line.strip() and line.replace(" ", "")[0] == "1":
            add_lines = True
        # Pretty much you add all the lines
        if add_lines:
            data_lines.append(line.strip())

    # Write the parsed data to the CSV file
    with open(output_csv_path, "a+") as fout:
        # Write the data rows
        for row in data_lines:
            print(row)
            fout.write(f"{protein_name}\t{ligand_name}\t{row}\n")


def convert_pdbqt_to_pdb(pdbqt_file: str, pdb_file: str, disable_bonding=False) -> None:
    """
    Convert a PDBQT file to a PDB file with Open Babel.

    :param pdbqt_file: path to the PDBQT input file
    :param pdb_file: path to the PDB output file
    :param disable_bonding: disable automatic bonding with Open Babel
    """
    cmd_args = [
        "obabel",
        "-ipdbqt",
        pdbqt_file,
        "-opdb",
        "-O",
        pdb_file,
    ]
    if disable_bonding:
        # "a" = read option
        # "b" = disable automatic bonding
        cmd_args += ["-ab"]

    cmd_return = subprocess.run(
        cmd_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    stdout = cmd_return.stdout.decode("utf-8")
    logging.debug(stdout)


real_number_pattern = r"[-+]?[0-9]*\.?[0-9]+(e[-+]?[0-9]+)?"
score_re = re.compile(rf"REMARK VINA RESULT:\s*(?P<affinity>{real_number_pattern})")


def merge_pdbqt(input_files: list, output_file_path: str):

    """
    Extract and save only lines starting with ATOM from multiple files,
    and combine them into one output file with a final TER line.

    Args:
        input_files (list): List of input file paths to be processed.
        output_file_path (str): Path to save the combined output file.
    """
    with open(output_file_path, "w") as outfile:

        for input_file_path in input_files:
            print(f"Merging {input_file_path}")
            with open(input_file_path, "r") as infile:
                for line in infile:
                    if line.startswith("ATOM"):
                        outfile.write(line)

        # Add the TER line to indicate termination
        outfile.write("TER\n")

    print(f"Combined ATOM lines saved to {output_file_path}")


def split_pdbqt(input_file, output_file_1, output_file_2):
    """
    Splits a PDBQT file into two parts based on the 'TER' separator and
    ensures the first output file always has more atoms than the second.

    Args:
        input_file (str): Path to the input PDBQT file.
        output_file_1 (str): Path to save the first split file.
        output_file_2 (str): Path to save the second split file.
    """
    with open(input_file, "r") as file:
        lines = file.readlines()

    # Split the lines based on 'TER'
    first_part = []
    second_part = []
    ter_encountered = False

    for line in lines:
        if line.strip() == "TER":
            if not ter_encountered:
                ter_encountered = True
                first_part.append(line)
            else:
                second_part.append(line)
        elif not ter_encountered:
            first_part.append(line)
        else:
            second_part.append(line)

    # Counting atoms in both parts
    num_atoms_1 = sum(1 for line in first_part if line.startswith("ATOM"))
    num_atoms_2 = sum(1 for line in second_part if line.startswith("ATOM"))

    # Ensure the first output file has more atoms
    if num_atoms_1 < num_atoms_2:
        first_part, second_part = second_part, first_part
        num_atoms_1, num_atoms_2 = num_atoms_2, num_atoms_1

    # Save both parts into separate files
    with open(output_file_1, "w") as file:
        file.writelines(first_part)

    with open(output_file_2, "w") as file:
        file.writelines(second_part)

    print(
        f"File split completed: {output_file_1} with {num_atoms_1} and {output_file_2} with {num_atoms_2} atoms."
    )

    return num_atoms_1, num_atoms_2


def get_coordinates_without_chain_id(pdb_file, seq_position):
    # Create a PDB parser
    parser = PDB.PDBParser(QUIET=True)

    # Parse the PDB file
    structure = parser.get_structure("protein", pdb_file)

    # Iterate over all chains and residues to find the matching sequence position
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[1] == seq_position:
                    print(residue)
                    # Extract the coordinates of the alpha carbon (CA) atom
                    ca_atom = residue["CA"]
                    return str(residue), ca_atom.get_coord()

    return None, None


def calculate_centroid(coords: list) -> tuple:

    """
    Calculate the geometric center (centroid) of a list of coordinates.

    Args:
        coords (list): List of XYZ coordinates.

    Returns:
        tuple: The XYZ coordinates of the centroid.
    """

    x_coords = [point[0] for point in coords]
    y_coords = [point[1] for point in coords]
    z_coords = [point[2] for point in coords]

    centroid = (
        sum(x_coords) / len(coords),
        sum(y_coords) / len(coords),
        sum(z_coords) / len(coords),
    )

    return centroid


####### Vina for apo ######## 


def mutate_and_save_pdb(parent_pdb, mutations, output_pdb):
    """
    Apply mutations to a PDB structure and save the mutated structure.

    Args:
        parent_pdb (str): Path to the parent PDB file.
        mutations (dict): Dictionary of mutations in the format {residue_id: new_aa}.
        output_pdb (str): Path to save the mutated PDB file.

    Returns:
        str: Path to the saved mutated PDB file.
    """

    # Load the parent structure
    universe = Universe(parent_pdb)

    # Apply all mutations
    for loc, aa in mutations.items():
        universe = apply_mutation(universe, (loc, aa))  # Ensure apply_mutation is defined

    # Save the mutated structure
    universe.atoms.write(output_pdb)
    print(f"Mutated structure saved to: {output_pdb}")

    return output_pdb



class VinaApoDock(ZSData):
    
    def __init__(
        self,
        input_csv: str,
        dock_opt: str,  #  ie "substrate", "joint", "all"
        cofactor_dets: str = "cofactor", # or inactivated_cofactor
        in_structure_dir: str = "data/structure",
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        fit_col_name: str = "fitness",
        output_dir: str = "zs/vina/apo",
        pH: float = 7.4,
        size_x: float = 20.0,
        size_y: float = 20.0,
        size_z: float = 20.0,
        num_modes: int = 9,
        exhaustiveness: int = 32,
        max_workers: int = 24,
        regen: bool = False,
        redock: bool = False
    ):

        super().__init__(
            input_csv=input_csv,
            combo_col_name=combo_col_name,
            fit_col_name=fit_col_name,
        )

        self._dock_opt = dock_opt
        self._cofactor_dets = cofactor_dets
        self._in_structure_dir = in_structure_dir
        self._output_dir = checkNgen_folder(os.path.join(output_dir, self.lib_name))
        self._common_pdbqt_dir = checkNgen_folder(os.path.join(self._output_dir, "common_pdbqt"))

        self._pH = pH
        self._size_x = size_x
        self._size_y = size_y
        self._size_z = size_z
        self._num_modes = num_modes
        self._exhaustiveness = exhaustiveness
        self._max_workers = max_workers

        self._regen = regen
        self._redock = redock

        self._coords = self._get_coords()

        self._ligand_pdbqt = os.path.join(self._common_pdbqt_dir, f"{self.substrate_dets}.pdbqt")
        self._cofactor2dock, self._cofactor2freeze = self._prep_common_pdbqt()

        self._dock_lib_parallel()


    def _get_coords(self):
        """
        Calculates the centroid of the ligand in the given structure.
        """

        coords = calculate_ligand_centroid(
            pdb_file=self.clean_struct, # the main pdb file dir
            ligand_info=ENZYME_INFO_DICT[self.protein_name]["ligand-info"]
        )

        return coords


    def _prep_common_pdbqt(self):
        """
        Prepares the PDBQT files for docking by converting smiles to PDBQT format.
        """

        # generate all pdbqt files from substrate and individual cofactor
        ligand_smiles2pdbqt(
            smiles=self.substrate_smiles,
            ligand_sdf_file=self._ligand_pdbqt.replace(".pdbqt", ".sdf"),
            ligand_pdbqt_file=self._ligand_pdbqt,
            pH=self._pH,
            regen=self._regen
        )

        cofactor2dock = []
        cofactor2freeze = []


        for (cofactor_name, cofactor_smiles) in zip(
            self.lib_info[self._cofactor_dets],
            self.lib_info[f"{self._cofactor_dets}-smiles"],
        ):  
            print(f"Processing cofactor {cofactor_name} {cofactor_smiles}")

            if self._dock_opt == "all":

                # check it is an ion
                simple_ion = cofactor_name.upper()[:2]
                if simple_ion in LIGAND_IONS or (cofactor_smiles and cofactor_smiles[0] == "[" and cofactor_smiles[-1] == "]"):

                    cofactor2dock.append(prepare_ion(
                        smiles=cofactor_smiles,
                        element=simple_ion if simple_ion in LIGAND_IONS else None,
                        output_dir=self._common_pdbqt_dir,
                        input_struct_path=self.clean_struct,
                        pH=self._pH,
                    ))

                else:
                    cofactor_pdbqt = os.path.join(self._common_pdbqt_dir, f"{cofactor_name}.pdbqt")
                    cofactor2dock.append(ligand_smiles2pdbqt(
                        smiles=cofactor_smiles,
                        ligand_sdf_file=cofactor_pdbqt.replace(".pdbqt", ".sdf"),
                        ligand_pdbqt_file=cofactor_pdbqt,
                        pH=self._pH,
                        regen=self._regen
                    ))
        
            # substrate only for trpb but sub + carbene for heme
            # extract heat atom to a pdb
            else:
                hetatm_pdb = os.path.join(self._common_pdbqt_dir, "hetatm.pdb")
                hetatm_pdbqt = os.path.join(self._common_pdbqt_dir, "hetatm.pdbqt")
                save_hetatm_only(
                    input_pdb=self.clean_struct,
                    output_pdb=hetatm_pdb
                )
                # obabel hetatm.pdb -O hetatm.pdbqt
                cmd = f"obabel {hetatm_pdb} -O {hetatm_pdbqt}"
                subprocess.run(cmd, shell=True)
                cofactor2freeze.append(hetatm_pdbqt)
            
        return cofactor2dock, cofactor2freeze
            

    def _mutate_apo(self, var):
        """
        Mutates the apo structure with the given mutation.
        """

        # make var_dir
        var_dir = checkNgen_folder(os.path.join(self._output_dir, var))
        var_mut_pdb = os.path.join(var_dir, f"{var}_mut.pdb")
        var_clean_pdb = os.path.join(var_dir, f"{var}.pdb")
        var_pdbqt = os.path.join(var_dir, f"{var}.pdbqt")

        if not self._regen and os.path.isfile(var_pdbqt):
            return var_pdbqt

        mutation_dict = {}

        if var != "WT":
            for v in var.split(":"):
                mutation_dict[int(v[1:-1])] = AA_DICT[v[-1]]

            # Mutate the apo structure
            mutate_and_save_pdb(
                parent_pdb=self.apo_struct,
                mutations=mutation_dict,
                output_pdb=var_mut_pdb
            )

        else:
            # copy the apo structure to the output directory
            shutil.copy(self.apo_struct, var_mut_pdb)
        
        # clean pdb and convert to pdbqt
        pdb_to_pdbqt_protein(var_mut_pdb, var_pdbqt)
        
        return var_pdbqt


    def _make_config(self, var):
        """
        Makes the config file for docking.
        """

        var_dir = checkNgen_folder(os.path.join(self._output_dir, var))
        conf_path = os.path.join(var_dir, f"{self._dock_opt}_conf.txt")

        # make the receptor pdbqt file
        if self._dock_opt != "all" and self._cofactor2freeze != []:
            # merge the pdbqt files
            receptor_pdbqt = os.path.join(var_dir, "receptor.pdbqt")
            merge_pdbqt(
                input_files=[self._mutate_apo(var)] + self._cofactor2freeze,
                output_file_path=receptor_pdbqt
            )

        else:
            receptor_pdbqt = self._mutate_apo(var)

        with open(conf_path, "w") as fout:
            fout.write(f"receptor = {receptor_pdbqt}\n")
            fout.write(f"ligand = {self._ligand_pdbqt}\n")

            # Include cofactors
            if self._cofactor2dock is not None:
                for cofactor_file in self._cofactor2dock:
                    fout.write(f"ligand = {cofactor_file}\n")

            fout.write(f"center_x = {self._coords[0]}\n")
            fout.write(f"center_y = {self._coords[1]}\n")
            fout.write(f"center_z = {self._coords[2]}\n")
            fout.write(f"size_x = {self._size_x}\n")
            fout.write(f"size_y = {self._size_y}\n")
            fout.write(f"size_z = {self._size_z}\n")
            fout.write(f"num_modes = {self._num_modes}\n")
            fout.write(f"exhaustiveness = {self._exhaustiveness}\n")

        return conf_path

    def _dock(self, var):
        """
        Dock the given variant.
        """

        # make the config file
        conf_path = self._make_config(var)
        docked_ligand_pdb = os.path.join(self._output_dir, var, f"{self._dock_opt}_docked.pdb")
        vina_logfile = os.path.join(self._output_dir, var, f"{self._dock_opt}_log.txt")

        # dock the variant
        cmd_list = ["vina", "--config", conf_path, "--out", docked_ligand_pdb, "--seed", str(42)]

        try:
            # Run the Vina command
            cmd_return = subprocess.run(
                cmd_list,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=True,  # Ensures subprocess raises an exception on error
            )

            # Write the Vina output to the log file
            with open(vina_logfile, "w") as fout:
                fout.write(cmd_return.stdout.decode("utf-8"))

        except subprocess.CalledProcessError as e:
            print(f"Error during Vina run: {e}")
            with open(vina_logfile, "w") as fout:
                fout.write(e.stdout.decode("utf-8"))
            raise
            
    def _dock_lib_parallel(self):
        """
        Dock the library in parallel.
        """

        # dock the library in parallel
        with ProcessPoolExecutor(max_workers=self._max_workers) as executor:
            futures = []
            for var in self.df[self._var_col_name]:
                if self._regen or self._redock:
                    futures.append(executor.submit(self._dock, var))
                else:
                    if not os.path.exists(os.path.join(self._output_dir, var, "log.txt")):
                        futures.append(executor.submit(self._dock, var))

            for future in as_completed(futures):
                future.result()


    @property
    def clean_struct(self):
        """The clean (no water or SO4) pdb file of the protein."""
        return os.path.join(self._in_structure_dir, "clean", f"{self.protein_name}.pdb")

    @property
    def apo_struct(self):
        """The apo pdb file of the protein."""
        return os.path.join(self._in_structure_dir, "apo", f"{self.protein_name}.pdb")
    
    @property
    def coords(self):
        return self._coords
            


###### extract docking scores ######


def extract_lowest_energy(file_path: str):

    with open(file_path, "r") as file:
        lines = file.readlines()

    for i, line in enumerate(lines):
        if "mode |   affinity | dist from best mode" in line:
            table_start = i + 3  # Skip to the actual table
            break

    energies = []
    for line in lines[table_start:]:
        if line.strip() == "":  # End of the table
            break
        parts = line.split()
        try:
            energies.append(float(parts[1]))  # Affinity is in the second column
        except ValueError:
            pass  # Skip invalid lines

    return min(energies)


def extract_binding_energy(log_file):
    """
    Extract the binding energy from an AutoDock Vina log file.

    Args:
        log_file (str): Path to the AutoDock Vina log file.

    Returns:
        float: Binding energy in kcal/mol, or None if not found.
    """
    with open(log_file, "r") as file:
        for line in file:
            if "Estimated Free Energy of Binding" in line:
                return float(line.split(":")[1].split()[0])


class VinaResults(ZSData):
    def __init__(
        self,
        input_csv: str,
        dock_opt: str,  #  ie "substrate",
        score_only: bool,  # = True,
        vina_struct_dir: str,  # = "vina/chai/struct_joint",
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        mut_col_name: str = "mut",
        pos_col_name: str = "pos",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        seq_dir: str = "data/seq",
        zs_dir: str = "zs",
        vina_score_dirname: str = "score",
        num_rep: int = 5,
        cofactor_type: str = "",
        withsub: bool = True,
    ):

        """
        A class to extract Vina docking energy.

        Args:
        - input_csv (str): The input CSV file, with `var_col_name` column
            ie. /disk2/fli/REVIVAL2/data/meta/not_scaled/PfTrpB-4bromo.csv
        - combo_col_name (str): The column name for the combo.
        - var_col_name (str): The column name for the variant.
        - mut_col_name (str): The column name for the mutation.
        - pos_col_name (str): The column name for the position.
        - seq_col_name (str): The column name for the sequence.
        - fit_col_name (str): The column name for the fitness.
        - seq_dir (str): The directory for the sequences.
        - zs_dir (str): The directory for the ZS data.
        - vina_dir (str): The directory for the Vina data.
        - vina_raw_dir (str): The directory for the raw Vina data.
        - vina_score_dir (str): The directory for the Vina score data.
        - num_rep (int): The number of replicates.
        - cofactor_type (str): The type of cofactor based on LIB_INFO_DICT.
            ie "" for simple cofacoctor, "inactivated" for inactivated cofactor, etc.
        - withsub (bool): Whether the zs includes info of the substrate.
        """

        super().__init__(
            input_csv=input_csv,
            combo_col_name=combo_col_name,
            var_col_name=var_col_name,
            mut_col_name=mut_col_name,
            pos_col_name=pos_col_name,
            seq_col_name=seq_col_name,
            fit_col_name=fit_col_name,
            withsub=withsub,
            seq_dir=seq_dir,
            zs_dir=zs_dir,
        )

        self._cofactor_append = f"_{cofactor_type}" if cofactor_type else ""
        self._vina_lib_name = f"{self.lib_name}{self._cofactor_append}"

        self._vina_dir = os.path.join(self._zs_dir, vina_struct_dir)
        self._vina_score_dir = checkNgen_folder(
            os.path.join(self._zs_dir, vina_struct_dir, vina_score_dirname)
        )
        self._vina_lib_dir = os.path.join(self._vina_dir, self._vina_lib_name)

        self._num_rep = num_rep
        self._dock_opt = dock_opt
        self._output_opt = "score_only" if score_only else "docked"
        print(f"Vina docking option: {self._dock_opt}")

        # for each row of the input csv, we will have a vina subfolder
        # ie I165A:I183A:Y301V_0 where 0 is the replicate number
        # with in the folder, there should be a file ends with _log.txt
        # which contains the docking energy
        # that can be extracted with extract_lowest_energy
        print(f"Extract vina score from {self._vina_lib_dir}")
        self._vina_df = self._extract_vina_score()

        # save the vina data
        print(f"Save vina score to {self.vina_df_path}")
        self._vina_df.to_csv(self.vina_df_path, index=False)


    def _extract_vina_score(self):

        """
        Extract the Vina data from the Vina log files.
        """

        df = self.df.copy()

        # Get the list of variants
        variants = df[self._var_col_name].unique()

        # Precompute scores
        scores = []
        for variant in tqdm(variants):
            for r in range(self._num_rep):
                vina_log_file = os.path.join(
                    self._vina_lib_dir,
                    f"{variant}_{r}",
                    f"{variant}_{r}-{self.substrate_name}-{self._dock_opt}-{self._output_opt}_log.txt",
                )

                # print(f"Extracting vina score for {variant}_{r} from {vina_log_file}")

                try:
                    if self._output_opt == "score_only":
                        vina_score = extract_binding_energy(vina_log_file)
                    else:
                        vina_score = extract_lowest_energy(vina_log_file)
                except Exception as e:
                    print(f"Error in extracting vina score for {variant}_{r}: {e}")
                    vina_score = np.nan

                scores.append((variant, r, vina_score))

        # Create a DataFrame from scores
        scores_df = pd.DataFrame(
            scores, columns=[self._var_col_name, "rep", "vina_score"]
        )

        # Pivot the scores DataFrame to match the original DataFrame structure
        scores_pivot = scores_df.pivot(
            index=self._var_col_name, columns="rep", values="vina_score"
        )
        scores_pivot.columns = [f"vina_{r}" for r in scores_pivot.columns]

        # Merge scores into the original DataFrame
        df = df.merge(scores_pivot, on=self._var_col_name, how="left")

        # take average of the replicates for each variant ignore nan when taking average
        df["vina"] = df[[f"vina_{r}" for r in range(self._num_rep)]].mean(axis=1)

        # add rank where the most negative score means rank 1
        df["vina_rank"] = df["vina"].rank(ascending=True)

        return df

    @property
    def vina_df(self):
        return self._vina_df

    @property
    def vina_df_path(self):
        return os.path.join(
            self._vina_score_dir,
            f"{self._vina_lib_name}-{self._dock_opt}-{self._output_opt}.csv",
        )

    @property
    def substrate_name(self):
        return LIB_INFO_DICT[self.lib_name]["substrate"]


def run_parse_vina_results(
    dock_opt: str,  #  ie "substrate",
    score_only: bool,  # = True,
    vina_struct_dir: str,  # = "vina/chai/struct_joint",
    pattern: Union[str, list] = "data/meta/not_scaled/*.csv",
    kwargs: dict = {},
):

    """
    Run the parse vina results function for all libraries

    Args:
    """
    if isinstance(pattern, str):
        lib_list = sorted(glob(pattern))
    else:
        lib_list = deepcopy(pattern)


    for lib in tqdm(lib_list):
        print(f"Running parse vina results for {lib}...")
        VinaResults(
            input_csv=lib,
            dock_opt=dock_opt,
            score_only=score_only,
            vina_struct_dir=vina_struct_dir,
            **kwargs,
        )