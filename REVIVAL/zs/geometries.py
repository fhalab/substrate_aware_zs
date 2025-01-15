"""
A script for calcalating the bond distances from generated structures
"""

from __future__ import annotations

import concurrent.futures

import os
from glob import glob
from tqdm import tqdm
from copy import deepcopy

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, MMCIFParser

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
    raise KeyError(
        f"Atom {atom_name} or its variations not found in residue {residue}."
    )


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
        raise ValueError(
            f"Unsupported geometry for atom with {len(neighbors)} neighbors."
        )

    print(
        f"Calculated hydrogen position from {atom} based on {neighbors}: {hydrogen_coord}"
    )

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
    add_hydrogen_to_2=False,
):
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

    if file_format.lower() == "pdb":
        parser = PDBParser(QUIET=True)
    elif file_format.lower() == "cif":
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

    print(
        f"Measuring bond distance between {chain_1}{residue_1}{atom_1} and {chain_2}{residue_2}{atom_2}."
    )

    # Determine coordinates for atom_1 and atom_2
    if add_hydrogen_to_1:
        neighbors_1 = get_covalent_neighbors(atom_1, residue_1)
        assert len(neighbors_1) <= 3, f"more than 3 neighbors: check{structure_file}"
        coord_1 = calculate_hydrogen_position(atom_1, neighbors_1)
    else:
        coord_1 = atom_1.coord

    if add_hydrogen_to_2:
        neighbors_2 = get_covalent_neighbors(atom_2, residue_2)
        assert len(neighbors_2) <= 3, f"more than 3 neighbors: check{structure_file}"
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
        struct_dir: str,
        scale_fit: str = "not_scaled",
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        mut_col_name: str = "mut",
        pos_col_name: str = "pos",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        seq_dir: str = "data/seq",
        zs_dir: str = "zs",
        bond_dir: str = "bonddist",
        max_workers: int = 64,
    ):

        """
        Initialize the BondData object.

        Args:
            input_csv (str): Path to the input CSV file.
            struct_dir (str): Path to the directory containing PDB or CIF files.
                ie. zs/af3/struct_joint/PfTrpB-4bromo
            scale_fit (str): Scale of the fitness values.
            combo_col_name (str): Column name for the combination of mutations.
            var_col_name (str): Column name for the variant name.
            mut_col_name (str): Column name for the mutation.
            pos_col_name (str): Column name for the position.
            seq_col_name (str): Column name for the sequence.
            fit_col_name (str): Column name for the fitness.
            seq_dir (str): Path to the directory containing sequence files.
            zs_dir (str): Path to the directory containing Z-score files.
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
            seq_dir=seq_dir,
            zs_dir=zs_dir,
        )

        self._max_workers = max_workers

        self._struct_dir = os.path.normpath(struct_dir)
        self._struct_subdir = os.path.join(struct_dir, self.lib_name)

        self._bond_dir = checkNgen_folder(
            os.path.join(
                zs_dir,
                bond_dir,
                os.path.basename(os.path.dirname(struct_dir)),
                os.path.basename(struct_dir),
            )
        )

        dock_opt = ""

        if "joint" in struct_dir:
            dock_opt = "_joint"
        elif "seperate" in struct_dir:
            dock_opt = "_seperate"

        if "chai" in struct_dir:
            dock_opt = "_chai" + dock_opt
        elif "af3" in struct_dir:
            dock_opt = "_af3" + dock_opt

        self._dist_opt = f"cofactor-distances{dock_opt}"

        # init the columns based on the keys from the dict
        # {'C-C': (('B', 1, 'LIG', 'C_5', False), ('B', 1, 'LIG', 'C_14', False)),
        # 'GLU-NH_1': (('A', 104, 'GLU', 'OE1', False), ('B', 1, 'LIG', 'N_1', True)),
        # 'GLU-NH_2': (('A', 104, 'GLU', 'OE2', False), ('B', 1, 'LIG', 'N_1', True))}
        self._dist_list = list(self.lib_info[self._dist_opt].keys())

        # both af and chai output 5 reps
        self._rep_list = [str(i) for i in list(range(5))] + ["avg"]

        print(f"self._dist_list: {self._dist_list}")
        print(f"self._rep_list: {self._rep_list}")

        print(f"Calculating bond distances for {self._struct_subdir}...")

        if "chai" in self._struct_subdir:
            df = self._append_chai_dist()

        elif "af3" in self._struct_subdir:
            df = self._append_af_dist()

        else:
            raise ValueError("Unsupported structure directory. Use 'chai' or 'af3'.")

        # save the dataframe
        output_csv = os.path.join(self._bond_dir, f"{self.lib_name}.csv")
        df.to_csv(output_csv, index=False)
        print(f"Saved bond distances to {output_csv}")

    def _calculate_distances(self, df, struct_type="af"):
        """
        Generalized method for calculating bond distances using ProcessPoolExecutor.
        """
        # Use external method for compatibility with ProcessPoolExecutor
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self._max_workers
        ) as executor:
            futures = [
                executor.submit(self._process_row, row, struct_type)
                for i, row in df.iterrows()
            ]

            # TQDM for progress tracking
            result_dfs = []
            for future in tqdm(
                concurrent.futures.as_completed(futures), total=len(futures)
            ):
                result_dfs.append(future.result())

        # Flatten results and convert to DataFrame
        flat_results = [item for sublist in result_dfs for item in sublist]
        results_df = pd.DataFrame(flat_results)

        # Pivot and merge
        df_wide = results_df.pivot_table(
            index=self._var_col_name, columns="key", values="distance", aggfunc="first"
        )

        return pd.merge(df, df_wide, on=self._var_col_name, how="outer")

    def _process_row(self, row, struct_type):

        """Helper method for parallel processing."""

        results = []
        if struct_type == "af":
            var_name = row[self._var_col_name].lower().replace(":", "_")
        else:
            var_name = row[self._var_col_name]

        var_struct_dir = os.path.join(self._struct_subdir, var_name)

        for d, dist in enumerate(self._dist_list):
            atom_info_1, atom_info_2 = self.lib_info[self._dist_opt][dist]

            atom_kwags = {
                "chain_id_1": atom_info_1[0],
                "res_id_1": atom_info_1[1],
                "atom_name_1": atom_info_1[3],
                "add_hydrogen_to_1": atom_info_1[4],
                "chain_id_2": atom_info_2[0],
                "res_id_2": atom_info_2[1],
                "atom_name_2": atom_info_2[3],
                "add_hydrogen_to_2": atom_info_2[4],
            }

            avg4agg = []

            for r in self._rep_list:

                # if not avg, then get the distance for each replicate
                if r != "avg":
                    structure_file = (
                        os.path.join(var_struct_dir, f"seed-1_sample-{r}", "model.cif")
                        if struct_type == "af"
                        else os.path.join(var_struct_dir, f"{var_name}_{r}.cif")
                    )

                    try:
                        distance = measure_bond_distance(
                            structure_file=structure_file, **atom_kwags
                        )
                        results.append(
                            {
                                self._var_col_name: row[self._var_col_name],
                                "key": f"{d}:{dist}_{r}",
                                "distance": distance,
                            }
                        )
                        avg4agg.append(distance)
                    except Exception as e:
                        print(f"Error processing {structure_file}: {e}")
                        results.append(
                            {
                                self._var_col_name: row[self._var_col_name],
                                "key": f"{d}:{dist}_{r}",
                                "distance": None,
                            }
                        )

                # now take the avg for both af and chai
                else:

                    avg_distance = np.mean(avg4agg) if avg4agg else None
                    results.append(
                        {
                            self._var_col_name: row[self._var_col_name],
                            "key": f"{d}:{dist}_avg",
                            "distance": avg_distance,
                        }
                    )

                    std_distance = np.std(avg4agg) if avg4agg else None
                    results.append(
                        {
                            self._var_col_name: row[self._var_col_name],
                            "key": f"{d}:{dist}_std",
                            "distance": std_distance,
                        }
                    )

            # add additional agg if af3
            if struct_type == "af":
                r = "agg"
                # ie zs/af3/struct_joint/PfTrpB-4bromo/i165a_i183a_y301v/i165a_i183a_y301v_model.cif
                structure_file = os.path.join(var_struct_dir, f"{var_name}_model.cif")
                try:
                    distance = measure_bond_distance(
                        structure_file=structure_file, **atom_kwags
                    )
                    results.append(
                        {
                            self._var_col_name: row[self._var_col_name],
                            "key": f"{d}:{dist}_{r}",
                            "distance": distance,
                        }
                    )
                    avg4agg.append(distance)
                except Exception as e:
                    print(f"Error processing {structure_file}: {e}")
                    results.append(
                        {
                            self._var_col_name: row[self._var_col_name],
                            "key": f"{d}:{dist}_{r}",
                            "distance": None,
                        }
                    )

        return results

    def _append_af_dist(self) -> pd.DataFrame:
        """Append distances for AF structures."""
        # df = self._init_dist_df()
        df = self.df.copy()
        return self._calculate_distances(df, struct_type="af")

    def _append_chai_dist(self) -> pd.DataFrame:
        """Append distances for CHAI structures."""
        # df = self._init_dist_df()
        df = self.df.copy()
        return self._calculate_distances(df, struct_type="chai")


def run_bonddist(
    struct_dir: str,
    pattern: str | list = "data/meta/not_scaled/*.csv",
    kwargs: dict = {},
):
    """
    Measure the bond distance for a library of structures.
    """

    if isinstance(pattern, str):
        lib_list = sorted(glob(pattern))
    else:
        lib_list = deepcopy(pattern)

    for lib in lib_list:
        print(f"Running parse vina results for {lib}...")
        BondData(input_csv=lib, struct_dir=struct_dir, **kwargs)