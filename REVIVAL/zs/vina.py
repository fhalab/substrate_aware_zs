"""
A script for get docking scores from chai structures

must use vina conda env
"""

from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Union

import logging
import subprocess

import warnings

import os
import re
from glob import glob
from tqdm import tqdm
from pathlib import Path
from copy import deepcopy

import numpy as np
import pandas as pd

from rdkit.Chem.MolStandardize.rdMolStandardize import Uncharger
from pdbfixer import PDBFixer
from openmm.app import PDBFile, PDBxFile
from Bio import PDB
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select, MMCIFIO

from rdkit import Chem
from rdkit.Chem import AllChem

from REVIVAL.global_param import LIB_INFO_DICT
from REVIVAL.preprocess import ZSData
from REVIVAL.util import checkNgen_folder, get_file_name, get_chain_structure, replace_residue_names_auto


warnings.filterwarnings("ignore")


def dock_lib_parallel(
    chai_dir: str,
    cofactor_type: str,
    residues: list = None,
    substrate_chain_ids: Union[list, str] = "B",
    freeze_opt: str = None,
    vina_dir: str = "vina",
    pH: float = 7.4,
    method="vina",
    size_x=10.0,
    size_y=10.0,
    size_z=10.0,
    num_modes=9,
    exhaustiveness=32,
    score_only=False,
    regen=False,
    rerun=False,
    max_workers=24,  # Number of parallel workers
):
    """
    A function to dock all generated chai structures and get scores in parallel.
    """

    output_dir = checkNgen_folder(chai_dir.replace("chai/mut_structure", vina_dir))

    keys_to_remove = [
        cofactor_type,
        cofactor_type.replace("-cofactor", ""),
        cofactor_type.replace("cofactor", ""),
        "cofactor",
        "-cofactor",
    ]
    lib_name = os.path.basename(chai_dir)
    for key in keys_to_remove:
        lib_name = lib_name.replace("_" + key, "")

    lib_dict = LIB_INFO_DICT[lib_name]

    cofactor_list = [
        (cofactor_smiles, cofactor, "B")
        for cofactor_smiles, cofactor in zip(
            lib_dict[cofactor_type + "-smiles"], lib_dict[cofactor_type]
        )
    ]

    print(f"Docking {chai_dir} with {cofactor_type} to {output_dir}")
    print(cofactor_list)
    print(f"Freeze option: {freeze_opt}")

    var_paths = sorted(glob(os.path.join(chai_dir, "*", "*.cif")))

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(
                dock_task,
                var_path=var_path,
                lib_dict=lib_dict,
                cofactor_list=cofactor_list,
                output_dir=output_dir,
                residues=residues,
                substrate_chain_ids=substrate_chain_ids,
                freeze_opt=freeze_opt,
                score_only=score_only,
                pH=pH,
                method=method,
                size_x=size_x,
                size_y=size_y,
                size_z=size_z,
                num_modes=num_modes,
                exhaustiveness=exhaustiveness,
                regen=regen,
                rerun=rerun,
            )
            for var_path in var_paths
        ]

        for future in tqdm(as_completed(futures), total=len(var_paths)):
            try:
                future.result()  # Raises an exception if the task failed
            except Exception as e:
                print(f"Task failed for {futures[future]}: {e}")


def dock_task(
    var_path,
    lib_dict,
    cofactor_list,
    output_dir,
    residues: list = None,
    substrate_chain_ids: Union[list, str] = "B",
    freeze_opt: str = None,
    score_only=False,
    pH=7.4,
    method="vina",
    size_x=10.0,
    size_y=10.0,
    size_z=10.0,
    num_modes=9,
    exhaustiveness=32,
    regen=False,
    rerun=False,
):
    """
    Task to dock a single .cif file.
    """

    if freeze_opt is None:
        append_name = ""
    else:
        append_name = freeze_opt

    log_txt = glob(
        os.path.join(output_dir, get_file_name(var_path), f"*{append_name}_log.txt")
    )

    var_dir = (
        os.path.normpath(
            checkNgen_folder(os.path.join(output_dir, get_file_name(var_path)))
        )
        + "/"
    )

    if len(log_txt) > 0 and (rerun is False):
        print(f"{var_path} already docked")
        print(f"See {log_txt}")
        return
    elif len(log_txt) == 0 and rerun is False:
        print(f"{var_path} not docked yet or successful")
        print("Rerun automatically set to True")

    try:
        dock(
            pdb_path=var_path,
            smiles=lib_dict["substrate-smiles"],
            ligand_name=lib_dict["substrate"],
            residues=residues,
            substrate_chain_ids=substrate_chain_ids,
            cofactors=cofactor_list,
            freeze_opt=freeze_opt,
            score_only=score_only,
            protein_dir=var_dir,
            ligand_dir=var_dir,
            output_dir=var_dir,
            pH=pH,
            method=method,
            size_x=size_x,
            size_y=size_y,
            size_z=size_z,
            num_modes=num_modes,
            exhaustiveness=exhaustiveness,
            regen=regen,
        )
    except Exception as e:
        print(f"Error in docking {var_path}: {e}")


def dock(
    pdb_path: str,
    smiles: str,
    ligand_name: str,
    cofactors: list,  # List of (smiles, name) tuples for cofactors
    protein_dir: str,
    ligand_dir: str,
    output_dir: str,
    pH: float,
    method: str,
    residues: list = None,
    substrate_chain_ids: Union[list, str] = "B",
    freeze_opt: str = None,
    score_only=False,
    size_x=10.0,
    size_y=10.0,
    size_z=10.0,
    num_modes=9,
    exhaustiveness=32,
    regen=False,
):
    """
    Dock a substrate with multiple cofactors to a protein using vina.
    """

    # Step 1: Process the protein
    protein_name, protein_pdb_file, protein_pdbqt = format_pdb(
        file_path=pdb_path, protein_dir=protein_dir, pH=pH, regen=regen
    )

    if score_only:
        from_pdb = True
        # TODO might need to extract the main ligand from the jointly docked file

    # Step 2: Process the main ligand
    ligand_pdbqt, _, ligand_sdf = format_ligand(
        smiles=smiles, name=ligand_name, ligand_dir=ligand_dir, pH=pH, regen=regen
    )

    # Step 3: Process cofactors
    cofactor_pdbqts = []
    for cofactor_smiles, cofactor_name, cofactor_chainid in cofactors:

        if freeze_opt in cofactor_name:
            from_pdb = True
        else:
            from_pdb = False

        cofactor_pdbqt, _, _ = format_ligand(
            smiles=cofactor_smiles,
            name=cofactor_name,
            ligand_dir=ligand_dir,
            pH=pH,
            full_pdb_path=pdb_path,
            chain_id=cofactor_chainid,
            regen=regen,
            from_pdb=from_pdb,
        )
        cofactor_pdbqts.append(cofactor_pdbqt)
    
    print(f"cofactor_pdbqt: {cofactor_pdbqt}")

    if method in ["vina", "ad4"]:
        dock_vina(
            input_pdb_path=pdb_path,
            ligand_pdbqt=ligand_pdbqt,
            protein_pdbqt=protein_pdbqt,
            ligand_name=ligand_name,
            protein_name=protein_name,
            output_dir=output_dir,
            residues=residues,
            substrate_chain_ids=substrate_chain_ids,
            cofactor_files=cofactor_pdbqts,
            freeze_opt=freeze_opt,
            size_x=size_x,
            size_y=size_y,
            size_z=size_z,
            num_modes=num_modes,
            exhaustiveness=exhaustiveness,
            method=method,
        )

    return None


def dock_vina(
    input_pdb_path: str,
    ligand_pdbqt: str,
    protein_pdbqt: str,
    ligand_name: str,
    protein_name: str,
    output_dir: str,
    residues: list = None,
    substrate_chain_ids: Union[list, str] = "B",
    cofactor_files: list = None,
    freeze_opt: str = None,
    size_x=10.0,
    size_y=10.0,
    size_z=10.0,
    num_modes=9,
    exhaustiveness=32,
    method="vina",
    num_cpus: int = None,
    seed=42,
):
    """Perform docking with Vina."""

    if freeze_opt is None:
        append_name = ""
    else:
        append_name = f"-{freeze_opt}"

    conf_path = os.path.join(
        output_dir, f"{protein_name}-{ligand_name}{append_name}_conf.txt"
    )

    # Step 1: Create config file
    make_config_for_vina(
        input_pdb_path=input_pdb_path,
        protein_pdbqt=protein_pdbqt,
        ligand_pdbqt=ligand_pdbqt,
        conf_path=conf_path,
        residues=residues,
        substrate_chain_ids=substrate_chain_ids,
        cofactor_files=cofactor_files,
        freeze_opt=freeze_opt,
        size_x=size_x,
        size_y=size_y,
        size_z=size_z,
        num_modes=num_modes,
        exhaustiveness=exhaustiveness,
    )

    # Auxiliary files
    vina_logfile = os.path.join(
        output_dir, f"{protein_name}-{ligand_name}{append_name}_log.txt"
    )
    docked_ligand_pdb = os.path.join(
        output_dir, f"{protein_name}-{ligand_name}{append_name}.pdb"
    )

    # Step 2: Perform docking
    cmd_list = [
        "vina",
        "--config",
        conf_path,
        "--out",
        docked_ligand_pdb,
        "--seed",
        str(seed),
    ]

    if num_cpus is not None:
        cmd_list += ["--cpu", str(num_cpus)]

    cmd_return = subprocess.run(
        cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    with open(vina_logfile, "w") as fout:
        fout.write(cmd_return.stdout.decode("utf-8"))

    # Step 3: Convert to PDB
    convert_pdbqt_to_pdb(
        pdbqt_file=docked_ligand_pdb, pdb_file=docked_ligand_pdb, disable_bonding=True
    )

    return vina_logfile


def make_config_for_vina(
    input_pdb_path: str,
    protein_pdbqt: str,
    ligand_pdbqt: str,
    conf_path: str,
    residues: list = None,
    substrate_chain_ids: Union[list, str] = "B",
    cofactor_files: list = None,
    freeze_opt: str = None,
    size_x=10.0,
    size_y=10.0,
    size_z=10.0,
    num_modes=9,
    exhaustiveness=32,
):
    """
    Create the config file for Vina, including cofactors if specified.

    Args:
    - input_pdb_path (str): Path to the input PDB file.
    - protein_pdbqt (str): Path to the protein PDBQT file.
    - ligand_pdbqt (str): Path to the ligand PDBQT file.
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

    main_dir = os.path.dirname(os.path.dirname(protein_pdbqt))

    if freeze_opt is None:
        receptor_pdbqt = protein_pdbqt
        cofactor2dock = deepcopy(cofactor_files)

    elif "cofactor" in freeze_opt.lower():
        # zs/vina/PfTrpB-4bromo_cofactor/WT_0/WT_0/WT_0.pdbqt will need to become
        # zs/vina/PfTrpB-4bromo_cofactor/WT_0/comb_cofactor.pdbqt
        receptor_pdbqt = os.path.join(main_dir, "comb_cofactor.pdbqt")
        receptor_pdb = os.path.join(main_dir, "comb_cofactor.pdb")

        # Combine the protein and cofactor files
        merge_pdbqt(
            input_files=[protein_pdbqt] + cofactor_files,
            output_file_path=receptor_pdbqt,
        )

        # Convert the combined PDBQT file to PDB
        convert_pdbqt_to_pdb(
            pdbqt_file=receptor_pdbqt, pdb_file=receptor_pdb, disable_bonding=True
        )

        cofactor2dock = None

    elif "bound-heme" in freeze_opt.lower():
        # just combine heme and Fe with ligand bound
        receptor_pdbqt = os.path.join(main_dir, "comb_bound-heme.pdbqt")
        receptor_pdb = os.path.join(main_dir, "comb_bound-heme.pdb")

        cofactor2freeze = []
        cofactor2dock = []

        for c in cofactor_files:
            keywords = ["bound-heme.pdbqt", "fe.pdbqt", "bound-heme", "bound heme", "heme-with", "heme with", "heme bound", "heme-bound"]
            if any(keyword in c.lower() for keyword in keywords):
                cofactor2freeze.append(c)
            else:
                cofactor2dock.append(c)

        print(f"cofactor2freeze: {cofactor2freeze}")
        print(f"cofactor2dock: {cofactor2dock}")

        # Combine the protein and cofactor files
        merge_pdbqt(
            input_files=[protein_pdbqt]+cofactor2freeze,
            output_file_path=receptor_pdbqt,
        )

        # Convert the combined PDBQT file to PDB
        convert_pdbqt_to_pdb(
            pdbqt_file=receptor_pdbqt, pdb_file=receptor_pdb, disable_bonding=True
        )

    elif "heme" in freeze_opt.lower():
        # just combine heme and Fe
        receptor_pdbqt = os.path.join(main_dir, "comb_heme.pdbqt")
        receptor_pdb = os.path.join(main_dir, "comb_heme.pdb")

        cofactor2freeze = []
        cofactor2dock = []

        for c in cofactor_files:
            keywords = ["bound-heme", "bound heme", "heme-with", "heme with", "heme bound", "heme-bound"]
            matched_keyword = next((keyword for keyword in keywords if keyword in c.lower()), None)
            
            if matched_keyword:
                heme_path = c.lower().replace(matched_keyword, "heme")
                heme_ligand_path = c.lower().replace(matched_keyword, "heme-ligand")

                checkNgen_folder(heme_path)
                checkNgen_folder(heme_ligand_path)

                # Need to split the heme and heme-ligand, always more atoms in the first file
                split_pdbqt(
                    input_file=c, 
                    output_file_1=heme_path, 
                    output_file_2=heme_ligand_path
                )
                cofactor2freeze.append(heme_path)
                
                matching_files = [c for c in cofactor_files if any(keyword in c.lower() for keyword in ["carbene"])]

                if len(matching_files) == 0:
                    cofactor2dock.append(heme_ligand_path)

            elif "fe.pdbqt" in c.lower():
                cofactor2freeze.append(c)
            else:
                cofactor2dock.append(c)

        print(f"cofactor2freeze: {cofactor2freeze}")
        print(f"cofactor2dock: {cofactor2dock}")

        # Combine the protein and cofactor files
        merge_pdbqt(
            input_files=[protein_pdbqt]+cofactor2freeze,
            output_file_path=receptor_pdbqt,
        )

        # Convert the combined PDBQT file to PDB
        convert_pdbqt_to_pdb(
            pdbqt_file=receptor_pdbqt, pdb_file=receptor_pdb, disable_bonding=True
        )

    else:
        raise ValueError("Invalid freeze option")

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


##### helper functions #####


def save_clean_protein(s, toFile, keep_chain="A", keep_all_protein_chains=True):
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


def remove_hydrogen_pdb(pdbFile, toFile):
    parser = (
        MMCIFParser(QUIET=True) if pdbFile[-4:] == ".cif" else PDBParser(QUIET=True)
    )
    s = parser.get_structure("x", pdbFile)

    class NoHydrogen(Select):
        def accept_atom(self, atom):
            if atom.element == "H" or atom.element == "D":
                return False
            return True

    io = MMCIFIO() if toFile[-4:] == ".cif" else PDBIO()
    io.set_structure(s)
    io.save(toFile, select=NoHydrogen())


def clean_one_pdb(proteinFile, toFile, keep_chain="keep_all"):
    """Clean and then remove stuff from a pdb file and then read-add in missing things."""
    if proteinFile[-4:] == ".cif":
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    s = parser.get_structure(proteinFile, proteinFile)
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

        # Now run
        clean_one_pdb(file_path, protein_pdb_file)

        # replace LIG_B etc with LIG
        replace_residue_names_auto(protein_pdb_file, protein_pdb_file)

        # Step 5: Convert to pdbqt for vina
        pdb_to_pdbqt_protein(
            input_path=protein_pdb_file, output_path=protein_pdbqt_file, pH=pH
        )

    return name, protein_pdb_file, protein_pdbqt_file


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


def format_ligand(
    smiles: str,
    name: str,
    ligand_dir: str,
    pH: float,
    from_pdb: bool = False,
    full_pdb_path=None,
    chain_id=None,
    regen=False,
):

    """
    Check if the ligand exists already; if not, handle formatting. Special handling for ions.
    """

    this_ligand_dir = checkNgen_folder(os.path.join(ligand_dir, name))

    ligand_pdbqt_file = os.path.join(this_ligand_dir, name + ".pdbqt")
    ligand_sdf_file = os.path.join(this_ligand_dir, name + ".sdf")

    fe_subdir = checkNgen_folder(os.path.join(ligand_dir, "Fe"))
    fe_pdb_file = os.path.join(fe_subdir, "Fe.pdb")
    fe_pdbqt_file = os.path.join(fe_subdir, "Fe.pdbqt")

    # Return existing files if already prepared
    if os.path.isfile(ligand_pdbqt_file) and not regen:
        return ligand_pdbqt_file, fe_pdbqt_file, ligand_sdf_file

    print(f"Processing {name}: {smiles}")

    try:
        # Handle ions separately
        if len(smiles) > 0 and ("[" == smiles[0]) and ("]" == smiles[-1]):  # Detect ions
            element_smiles = re.search(r"\[([A-Za-z]+)", smiles).group(1)
            print(f"Handling ion: {name} as {element_smiles}")
            print(f"extracting from {full_pdb_path}")
            # Generate PDB file for the ion
            ion_pdb_file = os.path.join(this_ligand_dir, name + ".pdb")
            extract_ions(
                input_file_path=full_pdb_path,
                output_file_path=ion_pdb_file,
                ion=element_smiles,
            )

            # Convert PDB to PDBQT using Open Babel
            subprocess.run(
                ["obabel", ion_pdb_file, "-O", ligand_pdbqt_file, "--metal"], check=True
            )

        else:
            # Standard ligand processing
            smiles = Chem.CanonSmiles(smiles, useChiral=True)
            protonated_smiles = protonate_smiles(smiles, pH=pH)
            protonated_mol = Chem.MolFromSmiles(protonated_smiles, sanitize=True)
            uncharger = Uncharger()
            mol = uncharger.uncharge(protonated_mol)
            mol = Chem.AddHs(mol)

            if mol.GetNumAtoms() == 0 or from_pdb:
                # TODO - Handle this case better
                print(f"Manually extract {smiles} of {name}")
                # Generate PDB file for the ion
                pdb_file = os.path.join(this_ligand_dir, name + ".pdb")
                print(f"from {full_pdb_path} to {pdb_file}")
                get_chain_structure(
                    input_file_path=full_pdb_path,
                    output_file_path=pdb_file,
                    chain_id=chain_id,
                )

                replace_residue_names_auto(pdb_file, pdb_file)

                processed_pdb_file = os.path.join(
                    this_ligand_dir, name + "_processed.pdb"
                )

                with open(pdb_file, "r") as infile, open(
                    fe_pdb_file, "w"
                ) as fe_out, open(processed_pdb_file, "w") as ligand_out:
                    for line in infile:
                        if line.startswith("HETATM") and "FE" in line:
                            fe_out.write(line)  # Write Fe atom to its PDB file
                        elif line.startswith("HETATM"):
                            # line = re.sub(r'LIG_C', 'LIG ', line)  # Replace "LIG_C" with "LIG"
                            ligand_out.write(
                                line
                            )  # Write other atoms to ligand PDB file
                    # Add TER and END for valid PDB files
                    fe_out.write("TER\nEND\n")
                    ligand_out.write("TER\nEND\n")

                if os.path.isfile(fe_pdb_file):
                    # Convert PDB to PDBQT using Open Babel
                    subprocess.run(
                        [
                            "obabel",
                            fe_pdb_file,
                            "-O",
                            fe_pdbqt_file,
                            "--metal"
                        ],
                        check=True,
                    )

                if from_pdb:
                        os.system(
                    f"obabel {processed_pdb_file} -xr -p {pH} --partialcharge gasteiger -O {ligand_pdbqt_file}"
                )
                else:

                    mol = Chem.MolFromPDBFile(processed_pdb_file, sanitize=True)
                    mol = Chem.AddHs(mol)
                    AllChem.EmbedMolecule(mol)
                    writer = Chem.SDWriter(ligand_sdf_file)
                    writer.write(mol)
                    writer.close()

                    os.system(
                        f"conda run -n vina mk_prepare_ligand.py -i {ligand_sdf_file} -o {ligand_pdbqt_file}"
                    )

            else:
                AllChem.EmbedMolecule(mol)
                writer = Chem.SDWriter(ligand_sdf_file)
                writer.write(mol)
                writer.close()

                os.system(
                    f"conda run -n vina mk_prepare_ligand.py -i {ligand_sdf_file} -o {ligand_pdbqt_file}"
                )

        # Validate the output
        if (
            not os.path.isfile(ligand_pdbqt_file)
            or os.path.getsize(ligand_pdbqt_file) == 0
        ):
            raise ValueError(f"PDBQT file for {name} is empty or missing.")
    except Exception as e:
        print(f"Error preparing ligand/ion: {name} - {str(e)}")
        raise

    return ligand_pdbqt_file, fe_pdbqt_file, ligand_sdf_file


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
                    if (
                        atom_name.upper() == ion.upper()
                        or residue_name.upper() == ion.upper()
                    ):
                        outfile.write(line)
                        ion_found = True

    if not ion_found:
        print(f"No ions matching '{ion}' were found in the file.")
    else:
        print(
            f"Ions matching '{ion}' were successfully extracted to {output_file_path}."
        )


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


def calculate_centroid(coords):
    x_coords = [point[0] for point in coords]
    y_coords = [point[1] for point in coords]
    z_coords = [point[2] for point in coords]

    centroid = (
        sum(x_coords) / len(coords),
        sum(y_coords) / len(coords),
        sum(z_coords) / len(coords),
    )

    return centroid


def calculate_chain_centroid(input_file, chain_ids):
    """
    Calculate the geometric center (centroid) of all atoms in the specified chain(s).

    Args:
        input_file (str): Path to the input PDB or CIF file.
        chain_ids (list of str): List of chain IDs to calculate the centroid for.

    Returns:
        tuple: The XYZ coordinates of the centroid.
    """
    # Determine the file type
    file_extension = os.path.splitext(input_file)[-1].lower()
    if file_extension == ".cif":
        parser = MMCIFParser(QUIET=True)
    elif file_extension == ".pdb":
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError(
            "Unsupported file format. Only PDB and CIF files are supported."
        )

    # Parse the structure
    structure = parser.get_structure("protein", input_file)

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
    with open(input_file, 'r') as file:
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
    with open(output_file_1, 'w') as file:
        file.writelines(first_part)

    with open(output_file_2, 'w') as file:
        file.writelines(second_part)

    print(f"File split completed: {output_file_1} with {num_atoms_1} and {output_file_2} with {num_atoms_2} atoms.")

    return num_atoms_1, num_atoms_2


###### extract docking scores ######


def extract_lowest_energy(file_path: str):
    with open(file_path, "r") as file:
        lines = file.readlines()

    # Initialize table_start outside the loop for safety
    table_start = None

    for i, line in enumerate(lines):
        if "mode |   affinity | dist from best mode" in line:
            table_start = i + 3  # Skip to the actual table
            break
    
    # If the header wasn't found, return np.nan and print a warning
    if table_start is None:
        print(f"Cannot find the table in {file_path}")
        return np.nan

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


class VinaResults(ZSData):
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
        vina_dir: str = "vina",
        vina_raw_dir: str = "",
        vina_score_dirname: str = "score",
        freeze_opt: str = None,
        num_rep: int = 5,
        cofactor_type: str = "",
        withsub: bool = True,
    ):

        """
        A class to extract Vina docking energy.

        Args:
        - input_csv (str): The input CSV file, with `var_col_name` column
            ie. /disk2/fli/REVIVAL2/data/meta/not_scaled/PfTrpB-4bromo.csv
        - scale_fit (str): The scaling of the fitness values.
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
        - vina_score_dirname (str): The directory for the Vina score data.
        - num_rep (int): The number of replicates.
        - cofactor_type (str): The type of cofactor based on LIB_INFO_DICT.
            ie "" for simple cofacoctor, "inactivated" for inactivated cofactor, etc.
        - withsub (bool): Whether the zs includes info of the substrate.
        """

        super().__init__(
            input_csv=input_csv,
            scale_fit=scale_fit,
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

        self._cofactor_append = f"_{cofactor_type}" if cofactor_type and cofactor_type != "cofactor" else ""
        
        self._freeze_opt_append = freeze_opt if not None else ""
        self._freeze_opt_subdir = f"freeze_{str(freeze_opt).lower()}"

        self._vina_lib_name = f"{self.lib_name}{self._cofactor_append}"

        self._vina_dir = os.path.join(self._zs_dir, vina_dir, vina_raw_dir)
        self._vina_score_dir = checkNgen_folder(
            os.path.join(self._zs_dir, vina_dir, vina_score_dirname, self._freeze_opt_subdir)
        )
        self._vina_lib_dir = os.path.join(self._vina_dir, self._vina_lib_name)

        self._num_rep = num_rep

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

        df = self.input_df.copy()

        # # init the vina score columns with nan
        # for r in range(self._num_rep):
        #     df[f"vina_{r}"] = np.nan

        # Get the list of variants
        variants = df[self._var_col_name].unique()

        # Precompute scores
        scores = []
        for variant in tqdm(variants):
            for r in range(self._num_rep):
                vina_log_files = glob(
                    os.path.join(self._vina_lib_dir, f"{variant}_{r}", f"*{self._freeze_opt_append}_log.txt")
                )

                if len(vina_log_files) >= 1 and os.path.isfile(vina_log_files[0]):
                    vina_score = extract_lowest_energy(vina_log_files[0])
                    if len(vina_log_files) > 1:
                        print(
                            f"Multiple log files found for {variant}_{r}, using the first one."
                        )
                else:
                    print(f"No log file found for {variant}_{r}")
                    vina_score = np.nan  # Default value if the file does not exist
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
        return os.path.join(self._vina_score_dir, f"{self._vina_lib_name}.csv")


def run_parse_vina_results(
    pattern: Union[str, list] = "data/meta/not_scaled/*.csv",
    cofactor_type: str = "",
    freeze_opt: str = None,
    kwargs: dict = {},
):

    """
    Run the parse vina results function for all libraries

    Args:
    """
    if isinstance(pattern, str):
        lib_list = glob(pattern)
    else:
        lib_list = deepcopy(pattern)

    print(lib_list)

    for lib in lib_list:
        print(f"Running parse vina results for {lib}...")
        VinaResults(
            input_csv=lib, 
            cofactor_type=cofactor_type,
            freeze_opt=freeze_opt,
            **kwargs
        )