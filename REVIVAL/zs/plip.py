# Modified from initial work by Lukas Radtke {lradtke@caltech.edu}
"""
Script for running PLIP on a PDB file.
"""

from __future__ import annotations

import os
import subprocess
from copy import deepcopy
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from glob import glob
from tqdm import tqdm

import numpy as np
import pandas as pd

from lxml import etree
import xml.etree.ElementTree as ET

from REVIVAL.preprocess import ZSData
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
    # if cif_or_pdb_file.endswith(".cif"):
    # obabel input.cif -O output.pdb clean up the cif file
    cmd = f"obabel {cif_or_pdb_file} -O {pdb_file} --remove HOH"
    subprocess.run(cmd, shell=True)

    # else:
    #     # Copy PDB directly
    #     checkNgen_folder(var_out_dir)
    #     os.system(f"cp {cif_or_pdb_file} {pdb_file}")

    # Run PLIP
    run_plip(pdb_file=pdb_file, output_dir=var_out_dir)


def run_lib_plip(
    in_dir: str, out_dir: str = "zs/plip", regen: bool = False, max_workers: int = 64
):

    """
    Get plip report for each of the variant in a given directory

    if  in_dir = 'data/structure'
        out_dir = 'zs/plip'
        will look for structure directly under the folder, i.e.
            data/structure/PfTrpB.pdb
            data/structure/Rma.pdb
        to generate plip results under pdb subdirectory in the out_dir, i.e.
            zs/plip/pdb/PfTrpB/
            zs/plip/pdb/Rma/

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
            var_out_dir = checkNgen_folder(
                os.path.join(out_dir, "pdb", "struct", variant_name)
            )
            tasks.append((file, var_out_dir, variant_name, regen))

    elif "af3" in in_dir:
        agg_cif_files = glob(f"{in_dir}/*/*_model.cif")
        rep_cif_files = glob(f"{in_dir}/*/*/model.cif")

        # Case 2: Nested folders with CIF files
        for cif_file in sorted(agg_cif_files + rep_cif_files):
            lib_name = os.path.basename(in_dir)
            struct_dets = in_dir.split("af3/")[-1].split(f"/{lib_name}")[0]

            lib_out_dir = checkNgen_folder(
                os.path.join(out_dir, "af3", struct_dets, lib_name)
            )
            variant_path = Path(cif_file).relative_to(Path(in_dir))
            variant_name = variant_path.parts[0].upper().replace("_", ":")

            if "_model.cif" in cif_file:
                rep_name = "agg"
            else:
                rep_name = variant_path.parts[1].split("sample-")[-1]

            var_out_dir = checkNgen_folder(
                os.path.join(lib_out_dir, f"{variant_name}_{rep_name}")
            )
            tasks.append((cif_file, var_out_dir, f"{variant_name}_{rep_name}", regen))

    elif "chai" in in_dir:
        # Case 3: Nested folders with CIF files
        for cif_file in sorted(glob(f"{in_dir}/**/*.cif")):
            lib_name = os.path.basename(in_dir)
            struct_dets = in_dir.split("chai/")[-1].split(f"/{lib_name}")[0]

            lib_out_dir = checkNgen_folder(
                os.path.join(out_dir, "chai", struct_dets, lib_name)
            )
            variant_name = get_file_name(cif_file)
            var_out_dir = checkNgen_folder(os.path.join(lib_out_dir, variant_name))
            tasks.append((cif_file, var_out_dir, variant_name, regen))

    # Parallelize the tasks using ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        list(tqdm(executor.map(process_task, tasks), total=len(tasks)))


def plip_xml2dict(xml_path: str) -> dict:

    """
    Parses an XML file and converts it to a dictionary.

    Args:
        xml_path (str): Path to the XML file.

    Returns:
        dict: Dictionary representation of the XML.
    """

    def xml2dict(element):
        """
        Recursively converts an lxml element to a dictionary.
        """
        # Base case: If the element has no children, return its text or attributes
        if not len(element) and not element.attrib:
            return element.text

        # Start the dictionary with attributes if present
        element_dict = {}
        if element.attrib:
            element_dict["@attributes"] = element.attrib

        # Process child elements
        for child in element:
            child_dict = xml2dict(child)
            if child.tag not in element_dict:
                # First occurrence of the tag
                element_dict[child.tag] = child_dict
            else:
                # Multiple occurrences of the same tag
                if not isinstance(element_dict[child.tag], list):
                    element_dict[child.tag] = [element_dict[child.tag]]
                element_dict[child.tag].append(child_dict)

        # Add element text if present and not just whitespace
        if element.text and element.text.strip():
            element_dict["#text"] = element.text.strip()

        return element_dict

    # Parse the XML file
    tree = etree.parse(xml_path)
    root = tree.getroot()

    # Convert the root element to a dictionary
    return xml2dict(root)


def get_plip_active_site_dict(xml_path: str) -> dict:
    """
    Extracts the active site residues from a PLIP XML file.

    Args:
        xml_path (str): Path to the PLIP XML file.

    Returns:
        dict: Dictionary of active site residues.
            ie {12: 'MET', 26: 'LEU', 37: 'GLY'}
    """

    # Parse the XML file
    plip_dict = plip_xml2dict(xml_path)

    # Extract the active site residues
    active_site_dict = {}

    bindingsite = plip_dict["bindingsite"]

    if isinstance(bindingsite, dict):
        for i in bindingsite["bs_residues"]["bs_residue"]:
            active_site_dict[int(i["#text"][:-1])] = i["@attributes"]["aa"]
    else:
        for i in bindingsite:
            for j in i["bs_residues"]["bs_residue"]:
                active_site_dict[int(j["#text"][:-1])] = j["@attributes"]["aa"]

    return active_site_dict


# Energy estimation functions
def estimate_hydrophobic_energy(distance):
    return -0.17 * (4 - distance) if distance < 4 else 0


def estimate_hydrogen_bond_energy(distance, donor_angle):
    return -1.5 * (2.5 - distance) * (donor_angle / 180) if distance < 2.5 else 0


def estimate_salt_bridge_energy(distance):
    return -2.0 * (4 - distance) if distance < 4 else 0


def estimate_metal_complex_energy(distance):
    return -5.0 * (2.5 - distance) if distance < 2.5 else 0


def estimate_water_bridge_energy(distance, donor_angle):
    return -1.0 * (2.8 - distance) * (donor_angle / 180) if distance < 2.8 else 0


def estimate_pi_stack_energy(distance):
    return -2.5 * (5 - distance) if distance < 5 else 0


def estimate_pi_cation_energy(distance):
    return -3.5 * (4 - distance) if distance < 4 else 0


def estimate_halogen_bond_energy(distance):
    return -1.5 * (3.5 - distance) if distance < 3.5 else 0


# Mapping dictionary
PLIP_INTERACTION_MAP = {
    "hydrophobic": ("hydrophobic_interaction", "dist", estimate_hydrophobic_energy),
    "hydrogen_bond": (
        "hydrogen_bond",
        "dist_d-a",
        estimate_hydrogen_bond_energy,
        "don_angle",
    ),
    "salt_bridge": ("salt_bridge", "dist", estimate_salt_bridge_energy),
    "metal_complex": ("metal_complex", "dist", estimate_metal_complex_energy),
    "water_bridge": (
        "water_bridge",
        "dist_d-a",
        estimate_water_bridge_energy,
        "don_angle",
    ),
    "pi_stack": ("pi_stack", "centdist", estimate_pi_stack_energy),
    "pi_cation": ("pi_cation_interaction", "dist", estimate_pi_cation_energy),
    "halogen_bond": ("halogen_bond", "dist", estimate_halogen_bond_energy),
}


# Parse PLIP XML Report
def parse_plip_report(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()

    interactions = {key: [] for key in PLIP_INTERACTION_MAP}

    for key, (tag, dist_key, func, *extra_keys) in PLIP_INTERACTION_MAP.items():
        for interaction in root.findall(f".//{tag}"):
            distance = float(interaction.find(dist_key).text)
            if extra_keys:
                extra_value = float(interaction.find(extra_keys[0]).text)
                interactions[key].append(func(distance, extra_value))
            else:
                interactions[key].append(func(distance))

    return interactions


# Calculate total stabilization energy
def calculate_total_energy(interactions):
    total_energy = sum(sum(values) for values in interactions.values())
    return total_energy


class PLIP_ZS(ZSData):
    """
    Class for running PLIP on a set of PDB files.
    """

    def __init__(
        self,
        input_csv: str,
        plip_dir: str,  # ie zs/plip/af3/struct_joint
        var_col_name: str = "var",
        fit_col_name: str = "fit",
    ):
        super().__init__(
            input_csv,
            var_col_name,
            fit_col_name,
        )

        self._plip_dir = plip_dir

        self._common_rep_list = [str(i) for i in range(5)]
        self._rep_list = (
            self._common_rep_list + ["agg"]
            if "af3" in plip_dir
            else self._common_rep_list
        )

        self._plip_rep_df = self._get_plip_rep_df()

        # save the df with rep
        self._plip_rep_df.to_csv(
            os.path.join(self.score_rep_dir, f"{self.lib_name}.csv"), index=False
        )

        self._plip_df = self._process_plip_rep_df()
        # save the df without rep
        self._plip_df.to_csv(
            os.path.join(self.score_dir, f"{self.lib_name}.csv"), index=False
        )

    def _get_plip_rep_df(self) -> pd.DataFrame:
        """
        Get the plip dataframe
        """

        # Get the list of plip xml files
        plip_xml_files = glob(
            f"{self._plip_dir}/{self.lib_name}/*/report.xml", recursive=True
        )

        # Create a dictionary to store the data
        df_list = []

        # Loop through each plip xml file and extract the data
        for xml_file in tqdm(plip_xml_files):

            # Get the variant name from the file path
            var_name, rep = xml_file.split("/")[-2].split("_")

            # Parse the XML file and extract interaction data
            interactions = parse_plip_report(xml_file)

            # Calculate total stabilization energy
            total_energy = calculate_total_energy(interactions)

            # Ensure all interaction keys exist in every row, initializing missing ones as empty lists
            interactions_complete = {
                key: list(interactions.get(key, [])) for key in PLIP_INTERACTION_MAP
            }

            # Append entry to df_list
            df_list.append(
                {
                    self._var_col_name: var_name,
                    "rep": rep,
                    "plip_naive_score": total_energy,
                    **interactions_complete,  # Ensures all expected keys are present
                }
            )

        # Convert the list of dictionaries to a pandas DataFrame
        plip_df = pd.DataFrame(df_list)

        print(plip_df.columns)

        # add number of interactions for each type and the sum of each type of interaction
        for key in PLIP_INTERACTION_MAP:
            plip_df[f"num_{key}"] = plip_df[key].apply(len)
            plip_df[f"sum_{key}"] = plip_df[key].apply(sum)

        # add total number of interactions
        plip_df["num_interactions"] = plip_df[
            [f"num_{key}" for key in PLIP_INTERACTION_MAP]
        ].sum(axis=1)

        return plip_df

    def _process_plip_rep_df(self) -> pd.DataFrame:

        """
        Process the plip dataframe with reps so that each row is for one variant
        and all the reps are appended after eaach column name
        Take the average of the reps for each variant from 0 to 4
        get the mean and std for each variant as additional columns
        if agg is one of the reps for each variant,
        append that for each variant but do not take that part of the avg
        """

        df = self._plip_rep_df.copy()
        cols_w_reps = [
            c
            for c in df.columns
            if c not in [self._var_col_name, "rep"] + list(PLIP_INTERACTION_MAP.keys())
        ]

        if "agg" in self._rep_list:
            # Separate rows where rep == "agg" for processing
            agg_rows = df[df["rep"] == "agg"].copy()

            # Append `_agg` to column names for agg rows
            agg_rows.rename(
                columns={col: f"{col}_agg" for col in cols_w_reps}, inplace=True
            )

        # Filter rows for rep in [0, 1, 2, 3, 4]
        filtered_df = df[df["rep"].isin(self._common_rep_list)]

        # Initialize a new DataFrame for results
        result = pd.DataFrame()

        # Process each column and compute values
        for col in cols_w_reps:
            # Pivot table to reshape data: `var` as index, `rep` values as columns
            reshaped = filtered_df.pivot(
                index=self._var_col_name, columns="rep", values=col
            )

            # Rename columns to include rep (e.g., pocket-plip-sasa_0, pocket-plip-sasa_1, ...)
            reshaped.columns = [f"{col}_{rep}" for rep in reshaped.columns]

            # Compute the average across rep values (ignoring NaN)
            reshaped[f"{col}_avg"] = reshaped.mean(axis=1)

            # Merge into the result DataFrame
            result = pd.concat([result, reshaped], axis=1)

        # Reset index for the result DataFrame
        result.reset_index(inplace=True)

        if "agg" in self._rep_list:
            # Merge back `agg_rows`
            merge_df = pd.merge(result, agg_rows, on=self._var_col_name, how="outer")
        else:
            merge_df = result

        # merge with the self.df to get fitness info
        return pd.merge(
            self.df[self.col2merge],
            merge_df,
            on=self._var_col_name,
            how="outer",
        )

    @property
    def col2merge(self) -> list:
        """
        Get the columns to merge
        """
        col2merge = [self._var_col_name, self._fit_col_name]

        if "selectivity" in self.df.columns:
            col2merge += ["selectivity"]

        return col2merge

    @property
    def score_dir(self) -> str:
        """
        Get the score directory
        """
        return checkNgen_folder(self._plip_dir.replace("struct", "score"))

    @property
    def score_rep_dir(self) -> str:
        """
        Get the score directory for rep
        """
        return checkNgen_folder(os.path.join(self.score_dir, "rep"))


def run_all_plip_zs(pattern: str | list, plip_dir: str, kwargs: dict = {}):

    if isinstance(pattern, str):
        lib_list = sorted(glob(pattern))
    else:
        lib_list = deepcopy(pattern)

    for lib in lib_list:
        PLIP_ZS(input_csv=lib, plip_dir=plip_dir, **kwargs)