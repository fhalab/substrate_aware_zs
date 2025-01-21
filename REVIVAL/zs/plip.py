"""Script for running PLIP on a PDB file."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

from tqdm import tqdm

from REVIVAL.util import checkNgen_folder, convert_cif_to_pdb


def run_plip(pdb_file, output_dir):
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


def run_lib_plip(in_dir: str, out_dir: str="zs/plip"):

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

    in_dir = Path(in_dir)
    out_dir = Path(out_dir)

    if in_dir.name == "structure":
        # Case 1: Directly under the folder
        for file in tqdm(in_dir.glob("*.pdb")):
            variant_name = file.stem
            lib_out_dir = out_dir / "holo" / variant_name

            # Use existing PDB file
            lib_out_dir.mkdir(parents=True, exist_ok=True)
            pdb_output_path = lib_out_dir / file.name
            pdb_output_path.write_bytes(file.read_bytes())

            # Run PLIP
            run_plip(pdb_output_path, lib_out_dir)

    elif "af3" in str(in_dir):
        # Case 2: Nested folders with CIF files
        for cif_file in tqdm(in_dir.glob("**/*.cif")):
            variant_path = cif_file.relative_to(in_dir)
            variant_name = variant_path.parts[0].upper().replace("_", ":")
            rep_name = variant_path.parts[-2:]  # Extract sample and seed information

            # Prepare output directory
            lib_out_dir = out_dir / "af3" / variant_path.parent
            if "agg" in rep_name:
                output_name = f"{variant_name}_agg"
            else:
                output_name = f"{variant_name}_{rep_name[-1].split('-')[-1]}"
            var_out_dir = lib_out_dir / output_name
            var_out_dir.mkdir(parents=True, exist_ok=True)

            # Convert CIF to PDB
            pdb_file = var_out_dir / cif_file.with_suffix(".pdb").name
            convert_cif_to_pdb(cif_file, str(pdb_file), ifsave=True)

            # Run PLIP
            run_plip(pdb_file, var_out_dir)

    elif "chai" in str(in_dir):
        # Case 3: Nested folders with CIF files
        for cif_file in tqdm(in_dir.glob("**/*.cif")):
            variant_path = cif_file.relative_to(in_dir)
            variant_name = variant_path.parts[0].upper().replace("_", ":")
            rep_name = variant_path.parts[-1].split("_")[-1]

            # Prepare output directory
            lib_out_dir = out_dir / "chai" / variant_path.parent
            output_name = f"{variant_name}_{rep_name}"
            var_out_dir = lib_out_dir / output_name
            var_out_dir.mkdir(parents=True, exist_ok=True)

            # Convert CIF to PDB
            pdb_file = var_out_dir / cif_file.with_suffix(".pdb").name
            convert_cif_to_pdb(cif_file, str(pdb_file), ifsave=True)

            # Run PLIP
            run_plip(pdb_file, var_out_dir)
    
    