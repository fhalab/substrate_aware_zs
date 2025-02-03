"""
Util functions for chem and MDanalysis
"""

from __future__ import annotations

import logging
import subprocess

from rdkit import Chem


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


def apply_mutation(universe, mutation):
    """
    Apply a mutation to a protein structure.

    Args:
        universe (MDAnalysis.Universe): The MDAnalysis Universe object.
        mutation (tuple): Mutation specified as (resid, new_resname).
                          Example: (56, "ALA")

    Returns:
        MDAnalysis.Universe: Updated Universe with the mutation applied.
    """
    resid, new_resname = mutation

    # Get the original residue name
    # original_resname = get_original_resname(universe, resid)

    # Select the residue to mutate
    residue = universe.select_atoms(f"resid {resid}").residues
    if len(residue) == 0:
        raise ValueError(f"Residue with resid {resid} not found in structure.")

    # Mutate the residue name
    residue.resnames = new_resname
    return universe
