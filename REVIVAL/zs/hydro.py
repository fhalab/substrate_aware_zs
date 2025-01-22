# File to automatically run hydrophobicity calculations as a ZS
# Author: Lukas Radtke {lradtke@caltech.edu}
# Date: 16/12/24
# Code reviewed by fzl

import numpy as np
import pandas as pd
import os 
import subprocess
import ast
import re
import xmltodict
from glob import glob
from copy import deepcopy
from tqdm import tqdm
# from scipy import stats
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors
from openbabel import pybel

from Bio.PDB import PDBParser

from REVIVAL.preprocess import ZSData
from REVIVAL.utils import AA_DICT, calculate_chain_centroid, get_protein_structure


HOPP_WOODS_SCALE = {
        'ALA': -0.500,
        'ARG':  3.000,
        'ASN':  0.200,
        'ASP':  3.000,
        'CYS': -1.000,
        'GLN':  0.200,
        'GLU':  3.000,
        'GLY':  0.000,
        'HIS': -0.500,
        'ILE': -1.800,
        'LEU': -1.800,
        'LYS':  3.000,
        'MET': -1.300,
        'PHE': -2.500,
        'PRO':  0.000,
        'SER':  0.300,
        'THR': -0.400,
        'TRP': -3.400,
        'TYR': -2.300,
        'VAL': -1.500,
    }

kyte_doolittle = {
    'ALA':  1.800,  
    'ARG': -4.500,  
    'ASN': -3.500,  
    'ASP': -3.500,  
    'CYS':  2.500,  
    'GLN': -3.500,  
    'GLU': -3.500,  
    'GLY': -0.400,  
    'HIS': -3.200,  
    'ILE':  4.500,  
    'LEU':  3.800,  
    'LYS': -3.900,  
    'MET':  1.900,  
    'PHE':  2.800,  
    'PRO': -1.600,  
    'SER': -0.800,  
    'THR': -0.700,  
    'TRP': -0.900,  
    'TYR': -1.300,  
    'VAL':  4.200,  
}

eisenberg_consensus = {
    'ALA':  0.620,
    'ARG': -2.530,
    'ASN': -0.780,
    'ASP': -0.900,
    'CYS':  0.290,
    'GLN': -0.850,
    'GLU': -0.740,
    'GLY':  0.480,
    'HIS': -0.400,
    'ILE':  1.380,
    'LEU':  1.060,
    'LYS': -1.500,
    'MET':  0.640,
    'PHE':  1.190,
    'PRO':  0.120,
    'SER': -0.180,
    'THR': -0.050,
    'TRP':  0.810,
    'TYR':  0.260,
    'VAL':  1.080,
}


def calc_hydrophobitcity(active_site, scale):
    """
    Returns the hydrophobicity score of a variant given the active site dictionary and a hydrophobicity scale. 
    """
    hydrophobicity = 0
    if isinstance(active_site, str):
        active_site = ast.literal_eval(active_site)
    for res in active_site[0].values():
        if len(res) != 1:
            hydrophobicity += scale[res]
        else:
            hydrophobicity += scale[AA_DICT[res]]
    hydrophobicity = hydrophobicity / len(active_site) 
    return hydrophobicity


def get_ligand_centroid(pdb_file: str, ligand_info: list):
    """
    Calculates the centroid of a given list of atoms specified by chain, residue, and atom names.

    Args:
        pdb_file (str): Path to the PDB file.
        ligand_info (list of tuples): List of atoms specified as (chain_id, residue_name, atom_name).

    Returns:
        tuple: Centroid coordinates as (x, y, z).
    """

    structure = get_protein_structure(pdb_file)

    atom_coords = []

    for chain_id, residue_name, atom_name in ligand_info:
        for model in structure:
            try:
                chain = model[chain_id]
                for residue in chain:
                    # Match residue name (flexible: partial match, case insensitive)
                    if residue_name.lower() in residue.resname.lower():
                        # Match atom name (flexible: ignore underscores and case differences)
                        for atom in residue:
                            if atom_name.replace("_", "").lower() == atom.name.replace("_", "").lower():
                                atom_coords.append(atom.coord)
            except KeyError:
                print(f"Chain {chain_id} not found in the structure.")
                continue

    if not atom_coords:
        raise ValueError("No matching atoms found in the structure.")

    # Calculate centroid
    return np.mean(atom_coords, axis=0)


def extract_active_site_by_radius(
    pdb_file: str,
    target_coord: np.array,
    target_chain="A",
    distance_threshold=10.0
    ):

    """
    Extracts a list of amino acids in the specified chain whose centroids (side chain or CA for glycine)
    are within a given distance from a specified (x, y, z) coordinate.

    Args:
        pdb_file (str): Path to the PDB file.
        target_coord (np.array): Target (x, y, z) coordinate.
        target_chain (str): Chain ID to search within (default is "A").
        distance_threshold (float): Distance threshold in Ångströms.

    Returns:
        list: A list of tuples containing residue information
              (e.g., [("GLY", 12), ("ALA", 25)]).
    """

    structure = get_protein_structure(pdb_file)

    nearby_residues = []

    # Iterate through all residues in the specified chain
    for model in structure:
        chain = model[target_chain]  # Access the specified chain
        for residue in chain:
            # Exclude backbone atoms (N, CA, C, O) and calculate centroid of side chain atoms
            side_chain_atoms = [atom for atom in residue if atom.name not in {"N", "CA", "C", "O"}]

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
                residue_info = (residue.resname, residue.id[1])
                nearby_residues.append(residue_info)

    return nearby_residues


class HydroData(ZSData):
    def __init__(
            self,
            input_csv = str,
            scale_fit: str = 'parent',
            combo_col_name: str = 'AAs',
            var_col_name: str = 'var',
            fit_col_name: str = 'fitness',
            seq_col_name: str = 'seq',
            structure_dir: str = 'data/structure',
            plip_dir: str = 'zs/plip',  
            active_site_radius: int = 10, 
        ):

        super().__init__(
            input_csv = input_csv,
            scale_fit = scale_fit,
            combo_col_name = combo_col_name,
            var_col_name = var_col_name,
            fit_col_name = fit_col_name,
            seq_col_name = seq_col_name,
            structure_dir = structure_dir,
        )

        self._plip_dir = plip_dir
        self._active_site_radius = active_site_radius

        # TODO

    # def extract_active_site_from_xml(self, xml_path):
    #     """
    #     Reads a .xml file generated by PLIP and extracts the active site residues.
    #     Return active site dictionary of format: 
    #     """
    #     with open(xml_path, 'r') as xml:
    #         report = xmltodict.parse(xml.read())
    #     bs_residues = report['report']['bindingsite']['bs_residues']['bs_residue']
    #     active_site_res_parent = {} #position: identity
    #     for bs_residue in bs_residues:
    #         if bs_residue['@contact'] == 'False':
    #             active_site_res_parent[bs_residue['#text'][:-1]] = bs_residue['@aa']
    #     return active_site_res_parent


    # def inference(self):
    #     # read the campaign data
    #     csv = pd.read_csv(self._input_csv, keep_default_na=False)
    #     csv = csv[~csv.apply(lambda row: row.astype(str).str.contains(r"\*").any(), axis=1)].reset_index(drop=True)  # filters for stop codons, marked as '*'
        
    #     # load parend pdb for campaign
    #     enzyme_name = csv[self.enzyme_col_name][1] # assumes its the same enzyme for every variant
    #     substrate_name = csv[self.substrate_col_name][1] # assumes its the same enzyme for every variant
    #     pdb = os.path.join(self.structure_dir, f'{enzyme_name}-{substrate_name}.pdb')

    #     # prepares mutated position for the radius-based active site definition
    #     campaign_name = os.path.basename(self._input_csv[:-4])

    #     if self.variant_structures_available == True:
    #         self.run_plip(csv)
        
    #     # initializing score_df with active sites 
    #     rows = []
    #     for _, variant in tqdm(csv.iterrows(), desc=f'Creating comparison df, reading in campaign {os.path.basename(self._input_csv)}...'):
            
    #         # generating radius-based active site 
    #         active_site_radius = self.extract_active_site_by_radius_parent(pdb, self.active_site_radius, self.mut_pos_list, variant)
            
    #         # generating plip-baed active site 
    #         if self.variant_structures_available == True:
    #             xml_path = os.path.join(self.results_dir, 'plip_reports', f'{variant[self.enzyme_col_name]}-{variant[self.substrate_col_name]}_{variant[self._var_col_name]}', 'report.xml')
    #             active_site_plip = self.extract_active_site_from_xml(xml_path)
    #             for i, mut_pos in enumerate(self.mut_pos_list):
    #                 active_site_plip[mut_pos] = variant[self._combo_col_name][i]  
    #         else:
    #             active_site_plip = {}

    #         gt_fitness = variant[self._fit_col_name]
    #         new_row = {'campaign_name': campaign_name, 'active_site_radius':[active_site_radius], 'active_site_plip':[active_site_plip], 'gt_fitness':gt_fitness, 'ZS1':[0.], 'ZS2':[0.]}
    #         rows.append(new_row)

    #     score_df = pd.DataFrame(data=rows)

    #     mol = self.extract_ligands_to_mol(pdb)

    #     # add the two ZS scores to the score_df
    #     for n, row in tqdm(score_df.iterrows(), desc=f'Calulating ZS for variants...'):
            
    #         # calculate hydrophobicity with hopp scale and 10A active site definition
    #         hydro_hopp_radius = calc_hydrophobitcity(row['active_site_radius'], HOPP_WOODS_SCALE)
    #         hydro_hopp_radius = (hydro_hopp_radius - (min(HOPP_WOODS_SCALE.values()))) / (max(HOPP_WOODS_SCALE.values()) - (min(HOPP_WOODS_SCALE.values())))
            
    #         # calculate hydrophobicity with hopp scale and PLIP active site definition
    #         if self.variant_structures_available == True:
    #             hydro_hopp_PLIP = calc_hydrophobitcity(row['active_site_plip'], HOPP_WOODS_SCALE)
    #         else:
    #             pass

    #         # calculate hydrophobicity (tpsa) of ligand
    #         tpsa = mol.calcdesc()["TPSA"]

    #         # calculate the ZS for radius-based active site definition
    #         ZS1 = tpsa - hydro_hopp_radius
    #         score_df.at[n, 'ZS1'] = ZS1

    #         # calculate the ZS for plip-based active site definition
    #         if self.variant_structures_available == True:
    #             ZS2 = tpsa - hydro_hopp_PLIP
    #         else:
    #             ZS2 = 0
    #         score_df.at[n, 'ZS2'] = ZS2

    #     # save the score df
    #     score_df_path = os.path.join(self.results_dir, f'{os.path.basename(self._input_csv)[:-4]}_score_df.csv')
    #     score_df.to_csv(score_df_path)
    #     print(f'ZS score .csv file saved to {score_df_path}.')

    @property
    def substrate_mol(self):
        return Chem.MolFromSmiles(self.lib_info["substrate-smiles"])

    @property
    def substrate_hydro_logp(self):
        return Descriptors.MolLogP(self.substrate_mol)

    @property
    def substrate_hydro_tpsa(self):
        return Chem.rdMolDescriptors.CalcTPSA(self.substrate_mol)

    @property
    def substrate_hydro_sasa(self):
        return Chem.MolSurf.LabuteASA(self.substrate_mol)



def run_hydro(pattern: str | list = None, kwargs: dict = {}):    
    
    if isinstance(pattern, str):
        lib_list = glob(pattern)
    else:
        lib_list = deepcopy(pattern)
    
    for lib in lib_list:
        HydroData(input_csv=lib, **kwargs)




