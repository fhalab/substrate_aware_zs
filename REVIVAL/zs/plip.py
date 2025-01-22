"""Script for running PLIP on a PDB file."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from glob import glob
from tqdm import tqdm

from REVIVAL.util import checkNgen_folder, get_file_name


def run_plip(pdb_file: str, output_dir: str):
    """
    Runs the PLIP command for a given PDB file and stores the results in the output directory.
    """

    checkNgen_folder(output_dir)

    # Define the log file path
    log_file = Path(output_dir) / "plip.log"

    cmd = [
        "python",
        "-m",
        "plip.plipcmd",
        "-f",
        os.path.abspath(pdb_file),
        "--out",
        os.path.abspath(output_dir),
        "--xml",
    ]
    
    # Run the command and redirect output to the log file
    try:
        with open(log_file, "w") as log:
            subprocess.run(cmd, check=True, stdout=log, stderr=log)
    except subprocess.CalledProcessError as e:
        print(f"PLIP execution failed for {pdb_file}. Check logs in {log_file}.")
        print(f"Error: {e}")

def process_task(task):
    """
    Processes a single task, which includes converting CIF to PDB and running PLIP.
    Args:
        task (tuple): Contains input file, output directory, and regen flag.
    """
    cif_or_pdb_file, var_out_dir, variant_name, regen = task

    var_xml = os.path.join(var_out_dir, "report.xml")

    # Check if output already exists and skip if regen is False
    if not regen and os.path.exists(var_xml):
        print(f"PLIP results for {var_xml} already exist. Skipping...")
        return

    # Prepare the PDB file path
    pdb_file = os.path.join(var_out_dir, f"{variant_name}.pdb")

    # Convert CIF to PDB if necessary
    if cif_or_pdb_file.endswith(".cif"):
        # obabel input.cif -O output.pdb
        cmd = f"obabel {cif_or_pdb_file} -O {pdb_file}"
        subprocess.run(cmd, shell=True)

    else:
        # Copy PDB directly
        checkNgen_folder(var_out_dir)
        os.system(f"cp {cif_or_pdb_file} {pdb_file}")

    # Run PLIP
    run_plip(pdb_file=pdb_file, output_dir=var_out_dir)


def run_lib_plip(in_dir: str, out_dir: str="zs/plip", regen: bool=False, max_workers: int = 64):

    """
    Get plip report for each of the variant in a given directory 

    if  in_dir = 'data/structure'
        out_dir = 'zs/plip'
        will look for structure directly under the folder, i.e.
            data/structure/PfTrpB.pdb
            data/structure/Rma.pdb
        to generate plip results under holo subdirectory in the out_dir, i.e.
            zs/plip/holo/PfTrpB/
            zs/plip/holo/Rma/

    if  in_dir = zs/af3/struct_joint/ParLQ
        out_dir = zs/plip
        will look for structures under the subfolders for each variant, i.e.
            zs/af3/struct_joint/ParLQ/w56e_y57k_l59f_q60d_f89w/w56e_y57k_l59f_q60d_f89w_model.cif
            zs/af3/struct_joint/ParLQ/w56e_y57k_l59f_q60d_f89w/seed-1_sample-0/model.cif
            zs/af3/struct_joint/ParLQ/w56e_y57k_l59f_q60d_f89w/seed-1_sample-1/model.cif
        to first convert cif to pdb and then
        to generate plip results under the out_dir that 
        perserve the structure details as well as consolidate and rename the variants and reps, i.e.
            zs/plip/af3/struct_joint/ParLQ/W56E:Y57K:L59F:Q60D:F89W_agg/
            zs/plip/af3/struct_joint/ParLQ/W56E:Y57K:L59F:Q60D:F89W_0/
            zs/plip/af3/struct_joint/ParLQ/W56E:Y57K:L59F:Q60D:F89W_1/

    if  in_dir = zs/chai/struct_joint/ParLQ
        out_dir = zs/plip
        will look for structures under the subfolders for each variant, i.e.
            zs/chai/struct_joint/ParLQ/W56A:Y57C:L59S:Q60E:F89G/W56A:Y57C:L59S:Q60E:F89G_0.cif
            zs/chai/struct_joint/ParLQ/W56A:Y57C:L59S:Q60E:F89G/W56A:Y57C:L59S:Q60E:F89G_1.cif
        to first convert cif to pdb and then
        to generate plip results under the out_dir that
        perserve the structure details, i.e.
            zs/plip/chai/struct_joint/ParLQ/W56A:Y57C:L59S:Q60E:F89G_0/
            zs/plip/chai/struct_joint/ParLQ/W56A:Y57C:L59S:Q60E:F89G_1/
    """

    in_dir = os.path.normpath(in_dir)
    out_dir = checkNgen_folder(out_dir)

    in_dir = os.path.normpath(in_dir)
    out_dir = checkNgen_folder(out_dir)

    tasks = []

    if os.path.basename(in_dir) == "structure":
        # Case 1: Directly under the folder
        for file in sorted(glob(f"{in_dir}/*.pdb")):
            variant_name = get_file_name(file)
            var_out_dir = checkNgen_folder(os.path.join(out_dir, "holo", variant_name))
            tasks.append((file, var_out_dir, variant_name, regen))

    elif "af3" in in_dir:
        agg_cif_files = glob(f"{in_dir}/*/*_model.cif")
        rep_cif_files = glob(f"{in_dir}/*/*/model.cif")

        # Case 2: Nested folders with CIF files
        for cif_file in sorted(agg_cif_files + rep_cif_files):
            lib_name = os.path.basename(in_dir)
            struct_dets = in_dir.split("af3/")[-1].split(f"/{lib_name}")[0]

            lib_out_dir = checkNgen_folder(os.path.join(out_dir, "af3", struct_dets, lib_name))
            variant_path = Path(cif_file).relative_to(Path(in_dir))
            variant_name = variant_path.parts[0].upper().replace("_", ":")

            if "_model.cif" in cif_file:
                rep_name = "agg"
            else:
                rep_name = variant_path.parts[1].split("sample-")[-1]

            var_out_dir = checkNgen_folder(os.path.join(lib_out_dir, f"{variant_name}_{rep_name}"))
            tasks.append((cif_file, var_out_dir, f"{variant_name}_{rep_name}", regen))

    elif "chai" in in_dir:
        # Case 3: Nested folders with CIF files
        for cif_file in sorted(glob(f"{in_dir}/**/*.cif")):
            lib_name = os.path.basename(in_dir)
            struct_dets = in_dir.split("chai/")[-1].split(f"/{lib_name}")[0]

            lib_out_dir = checkNgen_folder(os.path.join(out_dir, "chai", struct_dets, lib_name))
            variant_name = get_file_name(cif_file)
            var_out_dir = checkNgen_folder(os.path.join(lib_out_dir, variant_name))
            tasks.append((cif_file, var_out_dir, variant_name, regen))

    # Parallelize the tasks using ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        list(tqdm(executor.map(process_task, tasks), total=len(tasks)))