"""
A script for get docking scores from chai structures

must use vina conda env
"""

from concurrent.futures import ThreadPoolExecutor, as_completed

import logging
import subprocess

import warnings

import os
import re
from glob import glob
from tqdm import tqdm
from pathlib import Path

from rdkit.Chem.MolStandardize.rdMolStandardize import Uncharger
from pdbfixer import PDBFixer
from openmm.app import PDBFile, PDBxFile
from Bio import PDB
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select, MMCIFIO

from rdkit import Chem
from rdkit.Chem import AllChem

from REVIVAL.global_param import LIB_INFO_DICT
from REVIVAL.util import checkNgen_folder, get_file_name, get_chain_structure


warnings.filterwarnings("ignore")


def dock_task(
    var_path,
    lib_dict,
    cofactor_list,
    output_dir,
    pH,
    method,
    size_x,
    size_y,
    size_z,
    num_modes,
    exhaustiveness,
    rerun,
):
    """
    Task to dock a single .cif file.
    """
    log_txt = glob(os.path.join(output_dir, get_file_name(var_path), "*_log.txt"))
    if len(log_txt) > 0 and (rerun is False):
        print(f"{var_path} already docked")
        return

    var_dir = (
        os.path.normpath(
            checkNgen_folder(os.path.join(output_dir, get_file_name(var_path)))
        )
        + "/"
    )

    try:
        dock(
            pdb_path=var_path,
            smiles=lib_dict["substrate-smiles"],
            ligand_name=lib_dict["substrate"],
            residues=list(lib_dict["positions"].values()),
            cofactors=cofactor_list,
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
        )
    except Exception as e:
        print(f"Error in docking {var_path}: {e}")


def dock_lib_parallel(
    chai_dir: str,
    cofactor_type: str,
    vina_dir: str = "vina",
    pH: float = 7.4,
    method="vina",
    size_x=15.0,
    size_y=15.0,
    size_z=15.0,
    num_modes=9,
    exhaustiveness=32,
    rerun=False,
    max_workers=8,  # Number of parallel workers
):
    """
    A function to dock all generated chai structures and get scores in parallel.
    """
    output_dir = checkNgen_folder(chai_dir.replace("chai/mut_structure", vina_dir))
    lib_dict = LIB_INFO_DICT[os.path.basename(chai_dir).replace("-plp", "")]
    cofactor_list = [
        (cofactor_smiles, cofactor, "B")
        for cofactor_smiles, cofactor in zip(
            lib_dict[cofactor_type + "-smiles"], lib_dict[cofactor_type]
        )
    ]

    print(f"Docking {chai_dir} with {cofactor_type} to {output_dir}")
    print(cofactor_list)

    var_paths = sorted(glob(os.path.join(chai_dir, "*", "*.cif")))

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(
                dock_task,
                var_path=var_path,
                lib_dict=lib_dict,
                cofactor_list=cofactor_list,
                output_dir=output_dir,
                pH=pH,
                method=method,
                size_x=size_x,
                size_y=size_y,
                size_z=size_z,
                num_modes=num_modes,
                exhaustiveness=exhaustiveness,
                rerun=rerun,
            )
            for var_path in var_paths
        ]
        for future in tqdm(as_completed(futures), total=len(var_paths)):
            try:
                future.result()  # Raises an exception if the task failed
            except Exception as e:
                print(f"Task failed for {futures[future]}: {e}")


def dock_lib(
    chai_dir: str,
    cofactor_type: str,
    vina_dir: str = "vina",
    pH: float = 7.4,
    method="vina",
    size_x=15.0,
    size_y=15.0,
    size_z=15.0,
    num_modes=9,  # Dunno check vina docks using the defaut
    exhaustiveness=32,
    rerun=False,
):
    """
    A function for dock all generated chai structures and get scores

    Args:
    - chai_dir (str): The directory of chai structures.
        ie. /disk2/fli/REVIVAL2/zs/chai/mut_structure/PfTrpB-4bromo
            with the following structure:
                I165A:I183A:Y301V
                    I165A:I183A:Y301V_0.cif
                    I165A:I183A:Y301V_0.npz
                    ...
    - cofactor_type (str): The type of cofactor based on LIB_INFO_DICT
        ie. cofactor, inactivated-cofactor, etc.
    - vina_dir (str): The directory to save the vina results.
    - pH (float): The pH for docking.
    - method (str): The docking method, either 'vina' or 'ad4'.
    - size_x (float): The size in x for docking.
    - size_y (float): The size in y for docking.
    - size_z (float): The size in z for docking.
    - num_modes (int): The number of modes for docking.
    - exhaustiveness (int): The exhaustiveness for docking.
    - rerun (bool): Whether to rerun the docking.
    """

    output_dir = checkNgen_folder(chai_dir.replace("chai/mut_structure", vina_dir))

    lib_dict = LIB_INFO_DICT[os.path.basename(chai_dir).replace("-plp", "")]

    cofactor_list = []
    for cofactor_smiles, cofactor in zip(
        lib_dict[cofactor_type + "-smiles"], lib_dict[cofactor_type]
    ):
        cofactor_list.append((cofactor_smiles, cofactor, "B"))

    print(f"Docking {chai_dir} with {cofactor_type} to {output_dir}")
    print(cofactor_list)

    for var_path in tqdm(sorted(glob(os.path.join(chai_dir, "*", "*.cif")))):
        # ie zs/chai/mut_structure/PfTrpB-4bromo-plp/I165A:I183A:Y301V/I165A:I183A:Y301V_0.cif

        # check if the file is already docked
        log_txt = glob(os.path.join(output_dir, get_file_name(var_path), "*_log.txt"))
        if len(log_txt) > 0 and (rerun is False):
            print(f"{var_path} already docked")
            continue

        var_dir = (
            os.path.normpath(
                checkNgen_folder(os.path.join(output_dir, get_file_name(var_path)))
            )
            + "/"
        )

        try:
            dock(
                pdb_path=var_path,
                smiles=lib_dict["substrate-smiles"],
                ligand_name=lib_dict["substrate"],
                residues=list(lib_dict["positions"].values()),
                cofactors=cofactor_list,
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
            )
        except Exception as e:
            print(f"Error in docking {var_path}: {e}")


def dock(
    pdb_path: str,
    smiles: str,
    ligand_name: str,
    residues: list,
    cofactors: list,  # List of (smiles, name) tuples for cofactors
    protein_dir: str,
    ligand_dir: str,
    output_dir: str,
    pH: float,
    method: str,
    size_x=5.0,
    size_y=5.0,
    size_z=5.0,
    num_modes=9,
    exhaustiveness=32,
):
    """Dock a substrate with multiple cofactors to a protein using vina."""

    # Step 1: Process the protein
    protein_name, protein_pdb_file, protein_pdbqt = format_pdb(
        pdb_path, protein_dir, pH
    )

    # Step 2: Process the main ligand
    ligand_pdbqt, _, ligand_sdf = format_ligand(smiles, ligand_name, ligand_dir, pH)

    # Step 3: Process cofactors
    cofactor_pdbqts = []
    for cofactor_smiles, cofactor_name, cofactor_chainid in cofactors:
        print(f"Processing {cofactor_name}: {cofactor_smiles}")
        cofactor_pdbqt, _, _ = format_ligand(
            cofactor_smiles, cofactor_name, ligand_dir, pH, pdb_path, cofactor_chainid
        )
        cofactor_pdbqts.append(cofactor_pdbqt)

    if method in ["vina", "ad4"]:
        affinities = dock_vina(
            ligand_pdbqt=ligand_pdbqt,
            protein_pdbqt=protein_pdbqt,
            ligand_name=ligand_name,
            protein_name=protein_name,
            output_dir=output_dir,
            residues=residues,
            size_x=size_x,
            size_y=size_y,
            size_z=size_z,
            num_modes=num_modes,
            exhaustiveness=exhaustiveness,
            method=method,
            cofactor_files=cofactor_pdbqts,
        )

    return None


def dock_vina(
    ligand_pdbqt: str,
    protein_pdbqt: str,
    ligand_name: str,
    protein_name: str,
    output_dir: str,
    residues: list,
    size_x=5.0,
    size_y=5.0,
    size_z=5.0,
    num_modes=9,
    exhaustiveness=32,
    method="vina",
    num_cpus: int = None,
    seed=42,
    cofactor_files: list = None,
):
    """Perform docking with Vina."""
    conf_path = os.path.join(output_dir, f"{protein_name}-{ligand_name}_conf.txt")

    # Step 1: Create config file
    make_config_for_vina(
        protein_pdbqt,
        ligand_pdbqt,
        residues,
        conf_path,
        size_x=size_x,
        size_y=size_y,
        size_z=size_z,
        num_modes=num_modes,
        exhaustiveness=exhaustiveness,
        cofactor_files=cofactor_files,
    )

    # Auxiliary files
    vina_logfile = os.path.join(output_dir, f"{protein_name}-{ligand_name}_log.txt")
    docked_ligand_pdb = os.path.join(output_dir, f"{protein_name}-{ligand_name}.pdb")

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
    pdb_file: str,
    ligand_file: str,
    residues: list,
    output_file: str,
    size_x=5.0,
    size_y=5.0,
    size_z=5.0,
    num_modes=9,
    exhaustiveness=32,
    cofactor_files: list = None,
):
    """Create the config file for Vina, including cofactors if specified."""
    coords = []
    for position in residues:
        _, coordinates = get_coordinates_without_chain_id(pdb_file, int(position))
        if coordinates is not None:
            coords.append(coordinates)
    coords = calculate_centroid(coords)

    with open(output_file, "w") as fout:
        fout.write(f"receptor = {pdb_file}\n")
        fout.write(f"ligand = {ligand_file}\n")
        fout.write(f"center_x = {coords[0]}\n")
        fout.write(f"center_y = {coords[1]}\n")
        fout.write(f"center_z = {coords[2]}\n")
        fout.write(f"size_x = {size_x}\n")
        fout.write(f"size_y = {size_y}\n")
        fout.write(f"size_z = {size_z}\n")
        fout.write(f"num_modes = {num_modes}\n")
        fout.write(f"exhaustiveness = {exhaustiveness}\n")

        # Include cofactors
        if cofactor_files:
            for cofactor_file in cofactor_files:
                fout.write(f"ligand = {cofactor_file}\n")


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


def dock_ad4_pdbqt(
    ligand_pdbqt, protein_pdbqt, logfile, output_dir, protein_name, ligand_name
) -> None:
    """Run AD4.

    ../../software/x86_64Linux2/autogrid4 -p UYO78372.gpf -l UYO78372.glg
    pythonsh ../../enzymetk/prepare_dpf4.py -l ligand.pdbqt -r UYO78372.pdbqt -o UYO78372.dpf
    ../../software/x86_64Linux2/autodock4 -p UYO78372.dpf -l docking_results_UYO78372.dlg

    """
    package_root = Path(__file__).resolve().parent.parent

    gpf = os.path.join(output_dir, f"{protein_name}_{ligand_name}.gpf")
    glg = os.path.join(output_dir, f"{protein_name}_{ligand_name}.glg")
    dpf = os.path.join(output_dir, f"{protein_name}_{ligand_name}.dpf")
    dlg = os.path.join(output_dir, f"{protein_name}_{ligand_name}.dlg")

    os.chdir(output_dir)
    os.system(f'cp {protein_pdbqt} {os.path.join(output_dir, protein_name + ".pdbqt")}')
    os.system(f'cp {ligand_pdbqt} {os.path.join(output_dir, ligand_name + ".pdbqt")}')
    protein_pdbqt = protein_name + ".pdbqt"
    ligand_pdbqt = ligand_name + ".pdbqt"

    print(output_dir)
    # Step 1 prepare GPF
    cmd_list = [
        f"{package_root}/docko/deps/x86_64Linux2/mgltools_x86_64Linux2_1.5.7/bin/pythonsh",
        f"{package_root}/docko/deps/prepare_gpf.py",
        "-l",
        ligand_pdbqt,
        "-r",
        protein_pdbqt,
        "-o",
        gpf,
    ]
    print(" ".join(cmd_list))
    os.system(" ".join(cmd_list))

    # --------- Step 2 prepare GLG
    os.system(
        " ".join(
            [f"{package_root}/docko/deps/x86_64Linux2/autogrid4", "-p", gpf, "-l", glg]
        )
    )

    # --------- Step 3 prepare DPF
    cmd_list = [
        f"{package_root}/docko/deps/x86_64Linux2/mgltools_x86_64Linux2_1.5.7/bin/pythonsh",
        f"{package_root}/docko/deps/prepare_dpf4.py",
        "-l",
        ligand_pdbqt,
        "-r",
        protein_pdbqt,
        "-o",
        dpf,
    ]
    os.system(" ".join(cmd_list))

    # --------- FINALLY RUN AD4
    cmd_list = [
        f"{package_root}/docko/deps/x86_64Linux2/autodock4",
        "-p",
        dpf,
        "-l",
        dlg,
    ]

    os.system(" ".join(cmd_list))

    # They can get the results from here.
    return dlg


def dock_autodock_pdbqt(
    conf_path, log_path, out_path, seed, num_cpus: int = None
) -> None:
    """
    Run AutoDock Vina.

    :param conf_path: config
    :param log_path: path to log file
    :param out_path: path to output file
    :param seed: random seed
    :param num_cpus: number of CPU cores available to AutoDock Vina
    """
    cmd_list = [
        "vina",  # Needs to be installed as vina.
        "--config",
        conf_path,
        "--out",
        out_path,
        "--seed",
        str(seed),
    ]
    # ToDo: add in scoring function for ad4
    print(f"vina --config {conf_path} --out {out_path}")
    if num_cpus is not None:
        cmd_list += ["--cpu", str(num_cpus)]

    cmd_return = subprocess.run(
        cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    output = cmd_return.stdout.decode("utf-8")
    logging.debug(output)

    # Write output to the logging file
    with open(log_path, "w+") as fout:
        fout.write(output)

    # If failure, raise DockingError
    if cmd_return.returncode != 0:
        print(f"Docking with Vina failed: {output}")


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


def format_pdb(file_path: str, protein_dir: str, pH: float):

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

    # Now run
    clean_one_pdb(file_path, protein_pdb_file)

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
    full_pdb_path=None,
    chain_id=None,
):

    """
    Check if the ligand exists already; if not, handle formatting. Special handling for ions.
    """

    this_ligand_dir = checkNgen_folder(os.path.join(ligand_dir, name))

    ligand_pdbqt_file = os.path.join(this_ligand_dir, name + ".pdbqt")
    ligand_sdf_file = os.path.join(this_ligand_dir, name + ".sdf")

    fe_pdb_file = os.path.join(this_ligand_dir, "Fe.pdb")
    fe_pdbqt_file = os.path.join(this_ligand_dir, "Fe.pdbqt")

    # Return existing files if already prepared
    if os.path.isfile(ligand_pdbqt_file):
        return ligand_pdbqt_file, fe_pdbqt_file, ligand_sdf_file

    try:
        # Handle ions separately
        if ("[" == smiles[0]) and ("]" == smiles[-1]):  # Detect ions
            element_smiles = re.search(r"\[([A-Za-z]+)", smiles).group(1)
            print(f"Handling ion: {name} as {element_smiles}")
            # Generate PDB file for the ion
            ion_pdb_file = os.path.join(this_ligand_dir, name + ".pdb")
            with open(ion_pdb_file, "w") as pdb_file:
                pdb_file.write(
                    f"ATOM      1  {element_smiles}   ION     1       0.000   0.000   0.000  1.00  0.00      {element_smiles}\n"
                    "TER\n"
                    "END\n"
                )

            # Convert PDB to PDBQT using Open Babel
            subprocess.run(
                ["obabel", ion_pdb_file, "-O", ligand_pdbqt_file, "--metal"], check=True
            )
            # with open(ligand_pdbqt_file, "w") as fout:
            #     fout.write(f"ROOT\n")
            #     fout.write(f"ATOM      1  {element_smiles}  UNL     1       0.000   0.000   0.000  1.00  0.00      {element_smiles}\n")
            #     fout.write(f"ENDROOT\n")
            #     fout.write(f"TORSDOF 0\n")
        else:

            # Standard ligand processing
            smiles = Chem.CanonSmiles(smiles, useChiral=True)
            protonated_smiles = protonate_smiles(smiles, pH=pH)
            protonated_mol = Chem.MolFromSmiles(protonated_smiles, sanitize=True)
            uncharger = Uncharger()
            mol = uncharger.uncharge(protonated_mol)
            mol = Chem.AddHs(mol)

            if mol.GetNumAtoms() == 0:
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
                            "--partialcharge",
                            "gasteiger",
                        ],
                        check=True,
                    )

                mol = Chem.MolFromPDBFile(pdb_file, sanitize=True)
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