# File to automatically run hydrophobicity calculations as a ZS
# Modified by fzl from initial work by Lukas Radtke {lradtke@caltech.edu}
# Date: 16/12/24

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor

import warnings

import os
import math
from glob import glob
from copy import deepcopy
from tqdm import tqdm
import numpy as np
import pandas as pd

from scipy.spatial import ConvexHull


from rdkit import Chem
from rdkit.Chem import AllChem


import freesasa
from MDAnalysis import Universe
import MDAnalysis as mda
import tempfile

from substrate_aware.zs.plip import get_plip_active_site_dict
from substrate_aware.global_param import AA_DICT, ENZYME_INFO_DICT
from substrate_aware.preprocess import ZSData
from substrate_aware.chem_helper import smiles2mol
from substrate_aware.util import (
    calculate_chain_centroid,
    calculate_ligand_centroid,
    get_protein_structure,
    get_chain_ids,
    checkNgen_folder,
    get_file_name,
)


warnings.filterwarnings("ignore")


# Van der Waals radii for common atoms (in Å)
# VDW_RADII = {
#     "H": 1.20,  # Hydrogen
#     "C": 1.70,  # Carbon
#     "N": 1.55,  # Nitrogen
#     "O": 1.52,  # Oxygen
#     "S": 1.80,  # Sulfur
# }


# One-letter amino acid code to side chain SMILES
AA_SMILES = {
    "A": "C",  # Alanine (Methyl group)
    "R": "CCCNC(=N)N",  # Arginine (Guanidinium group)
    "N": "CC(=O)N",  # Asparagine (Amide group)
    "D": "CC(=O)O",  # Aspartate (Carboxylate group)
    "C": "CS",  # Cysteine (Thiol group)
    "E": "CCC(=O)O",  # Glutamate (Carboxylate group)
    "Q": "CCC(=O)N",  # Glutamine (Amide group)
    "G": "",  # Glycine (No side chain)
    "H": "CC1=CNC=N1",  # Histidine (Imidazole group)
    "I": "CCC(C)C",  # Isoleucine (Branched aliphatic)
    "L": "CC(C)CC",  # Leucine (Branched aliphatic)
    "K": "CCCCN",  # Lysine (Amine group)
    "M": "CCSC",  # Methionine (Thioether group)
    "F": "CC1=CC=CC=C1",  # Phenylalanine (Benzyl group)
    "P": "C1CCC1",  # Proline (Cyclic aliphatic)
    "S": "CO",  # Serine (Hydroxymethyl group)
    "T": "C(C)O",  # Threonine (Hydroxyethyl group)
    "W": "CC1=CNC2=CC=CC=C12",  # Tryptophan (Indole group)
    "Y": "CC1=CC=C(O)C=C1",  # Tyrosine (Phenol group)
    "V": "C(C)C",  # Valine (Isopropyl group)
}



def calculate_convex_hull_volume(mol):
    # Generate 3D conformer
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    conf = mol.GetConformer()
    
    # Extract atomic coordinates
    coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    
    # Compute the convex hull
    hull = ConvexHull(coords)
    return hull.volume  # Returns the volume of the convex hull



def calculate_smiles_vol(smiles: str) -> float:
    """
    Calculate the approximate volume of an amino acid side chain.

    Args:
        smiles (str): SMILES representation of the side chain.

    Returns:
        float: Volume of the side chain in Å³.
    """

    if smiles == "":
        return 0.0

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string.")

    # Add hydrogens
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())

    # Calculate the volume
    # volume = 0.0
    # for atom in mol.GetAtoms():
    #     symbol = atom.GetSymbol()
    #     radius = VDW_RADII.get(symbol, 1.50)  # Default radius for unknown atoms
    #     volume += (4 / 3) * math.pi * (radius ** 3)  # Volume of a sphere: (4/3)*π*r³

    return calculate_convex_hull_volume(mol)


AA_VOL = {
    k: calculate_smiles_vol(v) for k, v in AA_SMILES.items()
}


class VolData(ZSData):
    def __init__(
        self,
        input_csv: str,
        var_col_name: str = "var",
        combo_col_name: str = "AAs",
        fit_col_name: str = "fitness",
        subsmiles_col_name: str = "substrate-smiles",
        cofsmiles_col_name: str = "cofactor-smiles",
        active_site_radius: int = 10,
        vol_dir: str = "zs/vol",
        max_workers: int = 64,
        samesub: bool = True,
    ):

        super().__init__(
            input_csv=input_csv,
            var_col_name=var_col_name,
            combo_col_name=combo_col_name,
            fit_col_name=fit_col_name,
            subsmiles_col_name=subsmiles_col_name,
            cofsmiles_col_name=cofsmiles_col_name,
        )

        self._vol_dir = checkNgen_folder(vol_dir)

        print(f"Calculating naive volume for {self.lib_name}...")
        if samesub:
            vol_df = self._get_vol_df()
        else:
            vol_df = self._get_vol_df_diff_sub()

        merge_cols = [self._var_col_name, self._fit_col_name]

        if "selectivity" in self.df.columns:
            merge_cols.append("selectivity")
        # merge the volume DataFrame with the input DataFrame
        self._vol_df = pd.merge(
            self.df[merge_cols], vol_df, on=self._var_col_name, how="left"
        )

        # save the volume DataFrame to a CSV file
        self.vol_df.to_csv(self.vol_csv_path, index=False)

    def _calc_var_vol(self, var: str) -> float:
        """Calculate the naive volume of a variant"""

        # plus the volume of the parent AA
        # minus the volume of the variant AA

        var_vol = self.parent_vol

        if var == "WT":
            return var_vol

        for v in var.split(":"):

            parent_aa = v[0]
            mut_aa = v[-1]

            var_vol += AA_VOL[parent_aa] - AA_VOL[mut_aa]

        return var_vol

    def _get_vol_df(self) -> pd.DataFrame:
        """Get the naive volume DataFrame"""

        vol_df = []

        substrate_vol = calculate_smiles_vol(self.substrate_smiles)
        # cofactor_vol = calculate_smiles_vol(self.cofactor_smiles)
        # joint_vol = calculate_smiles_vol(self.joint_smiles)

        for var in tqdm(self.df[self._var_col_name]):
            var_vol = self._calc_var_vol(var)
            vol_df.append(
                {
                    self._var_col_name: var,
                    "var_vol": var_vol,
                    "substrate_vol": substrate_vol,
                    # "cofactor_vol": cofactor_vol,
                    # "joint_vol": joint_vol,
                    "var-substrate_vol": var_vol - substrate_vol,
                    # "var-cofactor_vol": var_vol - cofactor_vol,
                    # "var-joint_vol": var_vol - joint_vol,
                },
            )

        return pd.DataFrame(vol_df)

    def _get_vol_df_diff_sub(self) -> pd.DataFrame:
        """Get the naive volume DataFrame with different substrates"""

        vol_df = []

        for _, row in tqdm(self.df.iterrows()):
            var, substate_smiles, cofactor_smiles = (
                row[self._var_col_name], row[self._subsmiles_col_name], row[self._cofsmiles_col_name]
            )

            joint_smiles = substate_smiles + "." + cofactor_smiles
            var_vol = self._calc_var_vol(var)
            substrate_vol = calculate_smiles_vol(substate_smiles)
            # try:
            #     cofactor_vol = calculate_smiles_vol(cofactor_smiles)
            # except:
            #     cofactor_vol = np.nan
            # try:
            #     joint_vol = calculate_smiles_vol(joint_smiles)
            # except:
            #     joint_vol = np.nan

            vol_df.append(
                {
                    self._var_col_name: var,
                    "var_vol": var_vol,
                    "substrate_vol": substrate_vol,
                    # "cofactor_vol": cofactor_vol,
                    # "joint_vol": joint_vol,
                    "var-substrate_vol": var_vol - substrate_vol,
                    # "var-cofactor_vol": var_vol - cofactor_vol,
                    # "var-joint_vol": var_vol - joint_vol,
                }
            )

        return pd.DataFrame(vol_df)


    @property
    def parent_vol(self) -> list:
        """Get the parent active site volume"""
        return ENZYME_INFO_DICT[self.protein_name]["volume"]

    @property
    def substrate_smiles(self) -> str:
        """Get the substrate SMILES"""
        return self.lib_info["substrate-smiles"]

    @property
    def cofactor_smiles(self) -> str:
        """Get the cofactor SMILES"""
        return ".".join(self.lib_info["cofactor-smiles"])

    @property
    def joint_smiles(self) -> str:
        """Get the joint SMILES"""
        return self.substrate_smiles + "." + self.cofactor_smiles

    @property
    def substrate_vol(self) -> list:
        """Get the substrate volume"""
        return 

    @property
    def cofactor_vol(self) -> list:
        """Get the cofactor volume"""
        return calculate_smiles_vol(self.cofactor_smiles)

    @property
    def joint_vol(self) -> list:
        """Get the joint volume"""
        return calculate_smiles_vol(self.joint_smiles)

    @property
    def vol_df(self) -> pd.DataFrame:
        """Get the volume DataFrame"""
        return self._vol_df

    @property
    def vol_csv_path(self) -> str:
        """Get the volume CSV path"""
        return os.path.join(self._vol_dir, f"{self.lib_name}.csv")


def run_all_vol(pattern: str | list, samesub: bool = True, kwargs: dict = {}):

    if isinstance(pattern, str):
        lib_list = sorted(glob(pattern))
    else:
        lib_list = deepcopy(pattern)

    for lib in tqdm(lib_list):
        VolData(input_csv=lib, samesub=samesub, **kwargs)