# File to automatically run hydrophobicity calculations as a ZS
# Modified by fzl from initial work by Lukas Radtke {lradtke@caltech.edu}
# Date: 16/12/24

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor

import warnings

import os
from glob import glob
from copy import deepcopy
from tqdm import tqdm
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors, MolSurf

import freesasa
from MDAnalysis import Universe
import MDAnalysis as mda

from REVIVAL.zs.plip import get_plip_active_site_dict
from REVIVAL.preprocess import ZSData
from REVIVAL.util import (
    calculate_chain_centroid,
    calculate_ligand_centroid,
    get_protein_structure,
    get_chain_ids,
    smiles2mol,
    checkNgen_folder,
    get_file_name,
)


# warnings.filterwarnings("ignore")


KYTE_DOOLITTLE_SCALE = {
    "ALA": 1.800,
    "ARG": -4.500,
    "ASN": -3.500,
    "ASP": -3.500,
    "CYS": 2.500,
    "GLN": -3.500,
    "GLU": -3.500,
    "GLY": -0.400,
    "HIS": -3.200,
    "ILE": 4.500,
    "LEU": 3.800,
    "LYS": -3.900,
    "MET": 1.900,
    "PHE": 2.800,
    "PRO": -1.600,
    "SER": -0.800,
    "THR": -0.700,
    "TRP": -0.900,
    "TYR": -1.300,
    "VAL": 4.200,
}

HOPP_WOODS_SCALE = {
    "ALA": -0.500,
    "ARG": 3.000,
    "ASN": 0.200,
    "ASP": 3.000,
    "CYS": -1.000,
    "GLN": 0.200,
    "GLU": 3.000,
    "GLY": 0.000,
    "HIS": -0.500,
    "ILE": -1.800,
    "LEU": -1.800,
    "LYS": 3.000,
    "MET": -1.300,
    "PHE": -2.500,
    "PRO": 0.000,
    "SER": 0.300,
    "THR": -0.400,
    "TRP": -3.400,
    "TYR": -2.300,
    "VAL": -1.500,
}

EISENBERG_CONSENSUS_SCALE = {
    "ALA": 0.620,
    "ARG": -2.530,
    "ASN": -0.780,
    "ASP": -0.900,
    "CYS": 0.290,
    "GLN": -0.850,
    "GLU": -0.740,
    "GLY": 0.480,
    "HIS": -0.400,
    "ILE": 1.380,
    "LEU": 1.060,
    "LYS": -1.500,
    "MET": 0.640,
    "PHE": 1.190,
    "PRO": 0.120,
    "SER": -0.180,
    "THR": -0.050,
    "TRP": 0.810,
    "TYR": 0.260,
    "VAL": 1.080,
}


HYDRO_SCALES = {
    "kd": KYTE_DOOLITTLE_SCALE,
    "hw": HOPP_WOODS_SCALE,
    "ec": EISENBERG_CONSENSUS_SCALE,
}


ACTIVE_SITE_TYPES = ["pocket-plip", "pocket-subcentroid", "pocket-subcofcentroid"]


def calc_hydrophobitcity(active_site_dict: dict, hydro_scale: str) -> np.float:
    """
    Returns the hydrophobicity score of a variant given the active site dictionary and a hydrophobicity scale.

    Args:
        active_site_dict (dict): Dictionary of residues in the active site.
            ie. {12: 'MET', 26: 'LEU', 37: 'GLY'}
        hydro_scale (dict): Dictionary of amino acid hydrophobicity values.
            ie. kd

    Returns:
        float: Hydrophobicity score.
    """

    hydrophobicity = 0

    for (site, res) in active_site_dict.items():

        hydrophobicity += HYDRO_SCALES[hydro_scale][res]

    # normalize by the number of residues in the active site
    return hydrophobicity / len(active_site_dict)


def extract_active_site_by_radius(
    pdb_file: str, target_coord: np.array, target_chain="A", distance_threshold=10.0
) -> dict:

    """
    Extracts a list of amino acids in the specified chain whose centroids (side chain or CA for glycine)
    are within a given distance from a specified (x, y, z) coordinate.

    Args:
        pdb_file (str): Path to the PDB file.
        target_coord (np.array): Target (x, y, z) coordinate.
        target_chain (str): Chain ID to search within (default is "A").
        distance_threshold (float): Distance threshold in Ångströms.

    Returns:
        dict: A dictonary of tuples containing residue information
              (e.g., {12: "GLY", 25: "ALA"}).
    """

    structure = get_protein_structure(pdb_file)

    nearby_residues = {}

    # Iterate through all residues in the specified chain
    for model in structure:
        chain = model[target_chain]  # Access the specified chain
        for residue in chain:
            # Exclude backbone atoms (N, CA, C, O) and calculate centroid of side chain atoms
            side_chain_atoms = [
                atom for atom in residue if atom.name not in {"N", "CA", "C", "O"}
            ]

            if not side_chain_atoms:
                # Use the alpha carbon (CA) as the centroid for glycine
                if residue.resname == "GLY" and "CA" in residue:
                    ca_atom = residue["CA"]
                    centroid = np.array(ca_atom.coord)
                else:
                    # Skip residues with no side chains or CA
                    continue
            else:
                # Calculate the centroid of the side chain
                side_chain_coords = np.array([atom.coord for atom in side_chain_atoms])
                centroid = np.mean(side_chain_coords, axis=0)

            # Calculate distance between target coordinate and residue centroid
            distance = np.linalg.norm(centroid - target_coord)
            if distance <= distance_threshold:
                nearby_residues[int(residue.id[1])] = residue.resname

    return nearby_residues


def calculate_sasa_from_pdb(
    pdb_file, active_site_residues=None, chains=None, specific_atoms=None
):
    """
    Calculate SASA for active site residues, specified chains, or specific atoms.

    Args:
        pdb_file (str): Path to the PDB file containing the protein.
        active_site_residues (dict): Dict of active site residues as resid: resname.
        chains (list of str): List of chain identifiers to calculate SASA (e.g., ["B", "C"]).
        specific_atoms (list of tuple): List of specific atoms as (chain, residue, atom).

    Returns:
        float: SASA value for the specified selection.
    """
    # Ensure only one selection method is specified
    if sum(x is not None for x in [active_site_residues, chains, specific_atoms]) != 1:
        raise ValueError(
            "Specify only one of active_site_residues, chains, or specific_atoms."
        )

    # Load PDB file
    u = Universe(pdb_file)

    # Build selection query
    if active_site_residues:
        selection_query = " or ".join(
            [
                f"(resname {resname} and resid {resid} and segid A)"
                for resid, resname in active_site_residues.items()
            ]
        )
        append_temp = "activesite"
    elif chains:
        selection_query = " or ".join([f"segid {chain}" for chain in chains])
        append_temp = "chains"
    elif specific_atoms:
        # Generate relaxed, case-insensitive selection query for specific atoms
        broader_search_atoms = [
            f"(segid {chain.upper()} and resname {resname.upper()} and "
            f"(name {atom.upper()} or name {atom.replace('_', '').upper()} or "
            f"name {atom.replace('_', '').upper()[:-1]}_{atom[-1]}))"
            for chain, resname, atom in specific_atoms
        ]
        selection_query = " or ".join(broader_search_atoms)
        append_temp = "specificatoms"

    # Select atoms and validate
    selected_atoms = u.select_atoms(selection_query)
    if len(selected_atoms) == 0:
        raise ValueError("No atoms found for the specified selection.")

    # Write selected atoms to a temporary PDB file
    temp_file = f"{get_file_name(pdb_file)}_{append_temp}_temp_selection.pdb"
    with mda.Writer(temp_file, multiframe=False) as writer:
        writer.write(selected_atoms)

    # Calculate SASA
    structure = freesasa.Structure(temp_file)
    sasa = freesasa.calc(structure).totalArea()

    # Clean up
    os.remove(temp_file)

    return sasa


class HydroData(ZSData):
    def __init__(
        self,
        input_csv: str,
        plip_dir: str,  # ie zs/plip/af3/struct_joint
        scale_fit: str = "parent",
        var_col_name: str = "var",
        fit_col_name: str = "fitness",
        active_site_radius: int = 10,
        hydro_dir: str = "zs/hydro",
        max_workers: int = 64,
    ):

        super().__init__(
            input_csv=input_csv,
            scale_fit=scale_fit,
            var_col_name=var_col_name,
            fit_col_name=fit_col_name,
        )

        self._plip_dir = plip_dir
        self._hydro_dir = hydro_dir
        self._active_site_radius = active_site_radius
        self._max_workers = max_workers

        # TODO add opt for just pdb structure
        self._struct_dock_type = "af3" if "af3" in plip_dir else "chai"
        # rep list would be 0 to 4 for chai, 0 to 4 plus agg for af3
        self._common_rep_list = [str(i) for i in range(5)]
        self._rep_list = (
            self._common_rep_list + ["agg"]
            if "af3" in plip_dir
            else self._common_rep_list
        )

        self._struct_dock_dets = "joint" if "joint" in plip_dir else "separate"

        print(f"Calculating hydrophobicity for {self.lib_name}...")

        hydro_df = self._get_hydro_df()

        # proecess the hydro_df
        self._hydro_df = self._process_hydro_df(hydro_df)
        self._hydro_df.to_csv(self.hydro_df_path, index=False)

    def _process_hydro_df(self, df: pd.DataFrame) -> pd.DataFrame:

        if "agg" in self._rep_list:
            # Separate rows where rep == "agg" for processing
            agg_rows = df[df["rep"] == "agg"].copy()

            # Append `_agg` to column names for agg rows
            agg_rows.rename(
                columns={col: f"{col}_agg" for col in self.cols_w_reps}, inplace=True
            )

        # Filter rows for rep in [0, 1, 2, 3, 4]
        filtered_df = df[df["rep"].isin(self._common_rep_list)]

        # Initialize a new DataFrame for results
        result = pd.DataFrame()

        # Process each column and compute values
        for col in self.cols_w_reps:
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
            self.df[[self._var_col_name, self._fit_col_name]],
            merge_df,
            on=self._var_col_name,
            how="outer",
        )

    # go through each variant from self.df and replicate and calculate the hydrophobicity
    # make it parallelizable and then save the results to a dataframe
    def _get_hydro_df(self) -> pd.DataFrame:
        """
        Run the hydrophobicity calculations for each variant and replicate in parallel.
        """
        # Create a list of variant-replicate pairs
        var_rep_pairs = [
            (var, rep)
            for var in self.df[self._var_col_name].unique()
            for rep in self._rep_list
        ]

        # Use ProcessPoolExecutor for parallel processing
        with ProcessPoolExecutor(max_workers=self._max_workers) as executor:
            # Use executor.map with the worker function
            hydro_data = list(
                tqdm(
                    executor.map(self._hydro_worker, var_rep_pairs),
                    total=len(var_rep_pairs),
                    desc="Calculating hydrophobicity",
                )
            )

        # Separate successful and failed calculations
        successes = [entry for entry in hydro_data if "error" not in entry]
        failures = [entry for entry in hydro_data if "error" in entry]

        # Log failed pairs
        if failures:
            print(f"\nFailed calculations for {len(failures)} variant-replicate pairs:")
            for failure in failures:
                print(
                    f"Variant {failure['variant']}, Replicate {failure['replicate']}: {failure['error']}"
                )

        # Convert successful hydrophobicity data to a DataFrame
        return pd.DataFrame(successes)

    # def _get_hydro_df(self) -> pd.DataFrame:
    #     """
    #     Run the hydrophobicity calculations for each variant and replicate sequentially (no parallelism).
    #     This is useful for debugging.
    #     """
    #     # Create a list of variant-replicate pairs
    #     var_rep_pairs = [(var, rep) for var in self.df[self._var_col_name].unique() for rep in self._rep_list]

    #     # Initialize a list to store results
    #     hydro_data = []

    #     # Sequentially process each variant-replicate pair
    #     for var, rep in tqdm(var_rep_pairs, desc="Calculating hydrophobicity"):
    #         try:
    #             # Call the worker function directly for debugging
    #             result = self._hydro_worker((var, rep))
    #             hydro_data.append(result)
    #         except Exception as e:
    #             # Log the error
    #             print(f"Error processing Variant {var}, Replicate {rep}: {e}")
    #             hydro_data.append({"variant": var, "replicate": rep, "error": str(e)})

    #     # Convert the results into a DataFrame
    #     return pd.DataFrame(hydro_data)

    def _hydro_worker(self, pair):
        """
        Worker function for hydrophobicity calculation.
        """
        var, rep = pair
        try:
            return self._get_var_hydro(var, rep)
        except Exception as e:
            # Log the failure as a dictionary for easier debugging
            return {"variant": var, "replicate": rep, "error": str(e)}

    def _get_var_hydro(self, var: str, rep: str) -> dict:

        """
        Calculate the hydrophobicity of a variant
        """

        var_rep_name = f"{var}_{rep}"

        # zs/plip/af3/struct_joint/ParLQ/F89A_0/F89A_0.pdb
        var_path = os.path.join(
            self._plip_dir, self.lib_name, var_rep_name, var_rep_name + ".pdb"
        )
        # zs/plip/af3/struct_joint/ParLQ/F89A_0/report.xml
        xml_path = os.path.join(
            self._plip_dir, self.lib_name, var_rep_name, "report.xml"
        )

        hydro_dict = deepcopy(self.subcof_hydro_dict)

        if self._struct_dock_dets == "joint":
            ligand_chain_ids = ["B"]
            hydro_dict["substrate-sasa"] = calculate_sasa_from_pdb(
                pdb_file=var_path, specific_atoms=self.substrate_info
            )
            hydro_dict["substrate+cofactor-sasa"] = calculate_sasa_from_pdb(
                pdb_file=var_path, chains=ligand_chain_ids
            )
        else:
            ligand_chain_ids = get_chain_ids(var_path)[1:]
            hydro_dict["substrate-sasa"] = calculate_sasa_from_pdb(
                pdb_file=var_path, chains=[ligand_chain_ids[0]]
            )
            hydro_dict["substrate+cofactor-sasa"] = calculate_sasa_from_pdb(
                pdb_file=var_path, chains=ligand_chain_ids[1:]
            )

        # get the active site residues
        active_site_dict = {
            "pocket-plip": get_plip_active_site_dict(xml_path=xml_path),
            "pocket-subcentroid": extract_active_site_by_radius(
                pdb_file=var_path,
                target_coord=calculate_ligand_centroid(
                    pdb_file=var_path, ligand_info=self.substrate_info
                ),
                distance_threshold=self._active_site_radius,
            ),
            "pocket-subcofcentroid": extract_active_site_by_radius(
                pdb_file=var_path,
                target_coord=calculate_chain_centroid(
                    input_file=var_path, chain_ids=ligand_chain_ids
                ),
                distance_threshold=self._active_site_radius,
            ),
        }

        # first get the sasa
        for active_site_type, active_site_info in active_site_dict.items():
            hydro_dict[f"{active_site_type}-sasa"] = calculate_sasa_from_pdb(
                pdb_file=var_path, active_site_residues=active_site_info
            )

        # Calculate hydrophobicity with three different scale
        # nested for loop to calculate hydrophobicity for each active site method
        # and each hydrophobicity scale
        # key = f"active_site_dict[key]-{hydrophobicity scale option}"

        for key in active_site_dict.keys():
            for scale in HYDRO_SCALES.keys():
                hydro_dict[f"{key}-{scale}"] = calc_hydrophobitcity(
                    active_site_dict=active_site_dict[key], hydro_scale=scale
                )

        # also add in var and rep info
        hydro_dict[self._var_col_name] = var
        hydro_dict["rep"] = rep

        return hydro_dict

    @property
    def substrate_info(self) -> list:
        """Get the substrate atom from the library info"""
        return self.lib_info[f"{self.lib_info['substrate']}-info"]

    @property
    def substrate_mol(self) -> Chem.Mol:
        """Get the substrate molecule from the library info"""
        return smiles2mol(self.lib_info["substrate-smiles"])

    @property
    def subcof_jointsmiles(self) -> str:
        """Get the joint substrate-cofactor smiles from the library info"""
        return (
            self.lib_info["substrate-smiles"]
            + "."
            + ".".join(self.lib_info["cofactor-smiles"])
        )

    @property
    def subcof_jointmol(self) -> Chem.Mol:
        """Get the joint substrate-cofactor molecule from the library info"""
        return smiles2mol(self.subcof_jointsmiles)

    @property
    def subcof_hydro_dict(self) -> dict:
        """Get the hydrophobicity of the substrate and joint substrate-cofactor"""
        hydro_dict = {}
        for mol, mol_opt in zip(
            [self.substrate_mol, self.subcof_jointmol],
            ["substrate", "substrate+cofactor"],
        ):
            hydro_dict[f"{mol_opt}-logp"] = Descriptors.MolLogP(
                mol
            )  # Descriptors.MolLogP
            hydro_dict[f"{mol_opt}-tpsa"] = Descriptors.TPSA(mol)
            hydro_dict[f"{mol_opt}-asa"] = MolSurf.LabuteASA(mol)
        return hydro_dict

    @property
    def cols_w_reps(self) -> list:
        """Get the columns with replicate information"""
        return [
            f"{active_site_type}-{hydro}"
            for active_site_type in ACTIVE_SITE_TYPES
            for hydro in ["sasa"] + list(HYDRO_SCALES.keys())
        ]

    @property
    def hydro_df_path(self) -> str:
        """Get the path to the hydrophobicity DataFrame"""
        hydro_df_dir = checkNgen_folder(
            os.path.join(
                self._hydro_dir,
                self._struct_dock_type,
                f"struct_{self._struct_dock_dets}",
            )
        )
        return os.path.join(hydro_df_dir, f"{self.lib_name}.csv")

    @property
    def hydro_df(self) -> pd.DataFrame:
        """Get the hydrophobicity DataFrame"""
        return self._hydro_df


def run_all_hydro(pattern: str | list, plip_dir: str, kwargs: dict = {}):

    if isinstance(pattern, str):
        lib_list = sorted(glob(pattern))
    else:
        lib_list = deepcopy(pattern)

    for lib in lib_list:
        HydroData(input_csv=lib, plip_dir=plip_dir, **kwargs)