"""
A script for calcalating the bond distances from generated structures
"""

from __future__ import annotations

import os
import re
from tqdm import tqdm
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, MMCIFParser, PDBIO

from REVIVAL.preprocess import ZSData
from REVIVAL.util import checkNgen_folder, replace_residue_names_auto


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

        self._struct_dir = struct_dir
        self._bond_dir = checkNgen_folder(zs_dir, bond_dir)

        # init the columns based on the keys from the dict
        # {'C-C': (('B', 1, 'LIG', 'C_5', False), ('B', 1, 'LIG', 'C_14', False)),
        # 'GLU-NH_1': (('A', 104, 'GLU', 'OE1', False), ('B', 1, 'LIG', 'N_1', True)),
        # 'GLU-NH_2': (('A', 104, 'GLU', 'OE2', False), ('B', 1, 'LIG', 'N_1', True))}
        self._dist_list = list(self.lib_info["cofactor-distances"].keys())

        # both af and chai output 5 reps
        self._rep_list = [str(i) for i in list(range(5))] + ["agg"]

        if "chai" in struct_dir:
            output_csv = os.path.join(
                self._bond_dir, "chai", f"{self.lib_name}_bonddist.csv"
            )
            print(f"Calculating bond distances for {self._struct_dir}...")
            self.df = self._append_chai_dist()
            # save the dataframe
            self.df.to_csv(output_csv, index=False)
            print(f"Saved bond distances to {output_csv}")
        elif "af3" in struct_dir:
            output_csv = os.path.join(
                self._bond_dir, "af3", f"{self.lib_name}_bonddist.csv"
            )
            print(f"Calculating bond distances for {self._struct_dir}...")
            self.df = self._append_af_dist()
            # save the dataframe
            self.df.to_csv(output_csv, index=False)
            print(f"Saved bond distances to {output_csv}")
        else:
            raise ValueError("Unsupported structure directory. Use 'chai' or 'af3'.")

    def _init_dist_df(self) -> pd.DataFrame:

        """
        Initialize the dataframe with the columns from the input CSV file.
        """

        df = self.df.copy()

        # make a list of 0 to 4 in string and add "agg"
        # add the columns to the dataframe with nan values
        for i, dist in enumerate(self._dist_list):
            for j in self._rep_list:
                col_name = f"{i}:{dist}_{j}"
                df[col_name] = np.nan

        return df.copy()

    def _append_af_dist(self) -> pd.DataFrame:

        """
        Append the bond distances to the dataframe from AF structures.
        """

        df = self._init_dist_df()

        # loop over all variants and calculate the distances and dist_list and add them to the dataframe
        for i, row in tqdm(df.iterrows()):

            var_name = row[self.var_col_name].lower().replace(":", "_")

            # ie zs/af3/struct_joint/PfTrpB-4bromo/i165a_i183a_y301v
            struct_dir = os.path.join(self._struct_dir, var_name)

            for d, dist in enumerate(self._dist_list):
                # get the inputs
                atom_info_1, atom_info_2 = self.lib_info["cofactor-distances"][dist]
                chain_id_1, res_id_1, atom_name_1, add_hydrogen_to_1 = atom_info_1
                chain_id_2, res_id_2, atom_name_2, add_hydrogen_to_2 = atom_info_2

                for r in self._rep_list:
                    if r != "agg":
                        # ie zs/af3/struct_joint/PfTrpB-4bromo/i165a_i183a_y301v/seed-1_sample-0/model.cif
                        structure_file = os.path.join(
                            struct_dir, f"seed-1_sample-{r}", "model.cif"
                        )
                    else:
                        # zs/af3/struct_joint/PfTrpB-4bromo/i165a_i183a_y301v/i165a_i183a_y301v_model.cif
                        structure_file = os.path.join(
                            struct_dir, f"{var_name}_model.cif"
                        )

                    # get the distance
                    df.at[i, f"{d}:{dist}_{r}"] = measure_bond_distance(
                        structure_file=structure_file,
                        chain_id_1=chain_id_1,
                        res_id_1=res_id_1,
                        atom_name_1=atom_name_1,
                        chain_id_2=chain_id_2,
                        res_id_2=res_id_2,
                        atom_name_2=atom_name_2,
                        add_hydrogen_to_1=add_hydrogen_to_1,
                        add_hydrogen_to_2=add_hydrogen_to_2,
                    )

        return df

    def _append_chai_dist(self) -> pd.DataFrame:

        df = self._init_dist_df()

        # loop over all variants and calculate the distances and dist_list and add them to the dataframe
        for i, row in tqdm(df.iterrows()):

            var_name = row[self.var_col_name]

            # ie s/chai/mut_structure/PfTrpB-4bromo_cofactor/I165A:I183A:Y301V
            struct_dir = os.path.join(self._struct_dir, var_name)

            for d, dist in enumerate(self._dist_list):
                # get the inputs
                atom_info_1, atom_info_2 = self.lib_info["cofactor-distances"][dist]
                chain_id_1, res_id_1, atom_name_1, add_hydrogen_to_1 = atom_info_1
                chain_id_2, res_id_2, atom_name_2, add_hydrogen_to_2 = atom_info_2

                avg4agg = []

                for r in self._rep_list:
                    if r != "agg":
                        # ie zs/af3/struct_joint/PfTrpB-4bromo/i165a_i183a_y301v/seed-1_sample-0/model.cif
                        structure_file = os.path.join(struct_dir, f"{var_name}_{r}.cif")

                    # get the distance
                    bond_dist = measure_bond_distance(
                        structure_file=structure_file,
                        chain_id_1=chain_id_1,
                        res_id_1=res_id_1,
                        atom_name_1=atom_name_1,
                        chain_id_2=chain_id_2,
                        res_id_2=res_id_2,
                        atom_name_2=atom_name_2,
                        add_hydrogen_to_1=add_hydrogen_to_1,
                        add_hydrogen_to_2=add_hydrogen_to_2,
                    )

                    df.at[i, f"{d}:{dist}_{r}"] = bond_dist
                    avg4agg.append(bond_dist)

                if r == "agg":
                    df.at[i, f"{d}:{dist}_{r}"] = np.mean(avg4agg)

        return df


