import numpy as np
import pandas as pd
import subprocess
import os
import xmltodict
import ipdb
from tqdm import tqdm
import re
from rdkit import Chem
from rdkit.Chem import Descriptors
import math
from scipy import stats
import copy
from enzymeoracle.utils import dirpaths
import ast
import xml.etree.ElementTree as ET

from Bio.PDB import MMCIFParser, PDBIO

structure_dir = '/disk2/lukas/EnzymeOracle/data/struc_test/struc/chai_run_0'
campaign_dir = '/disk2/lukas/EnzymeOracle/data/struc_test/campaign/'
results_dir = '/disk2/lukas/EnzymeOracle/data/struc_test/results'

campaign_name = 'PfTrpB-5iodo'

aa_dict = {
    'A': 'ALA',  # Alanine
    'R': 'ARG',  # Arginine
    'N': 'ASN',  # Asparagine
    'D': 'ASP',  # Aspartic acid
    'C': 'CYS',  # Cysteine
    'E': 'GLU',  # Glutamic acid
    'Q': 'GLN',  # Glutamine
    'G': 'GLY',  # Glycine
    'H': 'HIS',  # Histidine
    'I': 'ILE',  # Isoleucine
    'L': 'LEU',  # Leucine
    'K': 'LYS',  # Lysine
    'M': 'MET',  # Methionine
    'F': 'PHE',  # Phenylalanine
    'P': 'PRO',  # Proline
    'S': 'SER',  # Serine
    'T': 'THR',  # Threonine
    'W': 'TRP',  # Tryptophan
    'Y': 'TYR',  # Tyrosine
    'V': 'VAL',  # Valine
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

def extract_ligands_to_mol(pdb_file):
    hetatm_lines = []
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith('ATOM') and 'LIG' in line or line.startswith('HETATM'):
                hetatm_lines.append(line)
    with open('temp_hetatm.pdb', 'w') as temp_file:
        temp_file.writelines(hetatm_lines)
    mol = Chem.MolFromPDBFile('temp_hetatm.pdb')
    return mol

def extract_mutated_residues(seq_csv, row=1):
    mutation_str = seq_csv.iloc[row, 2]
    mutated_pos = re.findall(r'\d+', mutation_str)
    mutated_pos = ['A' + num for num in mutated_pos]
    mutated_pos = " ".join(mutated_pos)
    return mutated_pos


def extract_active_site_from_xml(xml_path):
    with open(xml_path, 'r') as xml:
        report = xmltodict.parse(xml.read())
    bs_residues = report['report']['bindingsite']['bs_residues']['bs_residue']
    active_site_res_parent = {} #position: identity
    for bs_residue in bs_residues:
        if bs_residue['@contact'] == 'False':
            active_site_res_parent[bs_residue['#text'][:-1]] = bs_residue['@aa']
    return active_site_res_parent


def extract_ca_coordinates(pdb_file):
    coordinates = {}
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith('ATOM'):
                atom_name = line[12:16].strip()
                res_index = int(line[22:26].strip())
                if atom_name == 'CA':
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coordinates[res_index] = [x, y, z]
    return coordinates


def extract_hetatm_coordinates(pdb_file):
    coordinates = []
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith('HETATM') or 'LIG' in line:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coordinates.append([x, y, z])
    return coordinates

def obtain_residue_from_index(pdb_file, idx):
    with open(pdb_file, 'r') as PDB:
        for line in PDB:
            if line.startswith("ATOM"): 
                res_seq = line[22:26].strip()  
                atom_name = line[12:16].strip()  
                
                if res_seq == str(idx) and atom_name == 'CA':
                    residue = line[17:20].strip()  
                    return residue
    return None  

def extract_active_site_by_radius(pdb_path, radius_A):
    hetatm_coords = extract_hetatm_coordinates(pdb_path)
    ca_coords_idx = extract_ca_coordinates(pdb_path)
    ca_coords = np.array(list(ca_coords_idx.values()))
    ca_indices = np.array(list(ca_coords_idx.keys()))
    active_site_ca_idx = []
    for hetatm_coord in hetatm_coords:
        hetatm_coords = np.array(hetatm_coord)
        indices_within_radius = np.where(np.linalg.norm(ca_coords - hetatm_coord, axis=1) <= radius_A)[0]
        active_site_ca_idx.extend(ca_indices[indices_within_radius])
    active_site_ca_idx = set(active_site_ca_idx)
    active_site = {}
    for active_site_idx in active_site_ca_idx:
        active_site[str(active_site_idx)] = obtain_residue_from_index(pdb_path, active_site_idx)
    return active_site


def is_number(s):
    try:
        float(s)
        return True
    except (ValueError, TypeError):
        return False


run_plip_flag = False
if run_plip_flag == True:
    for campaign_name in os.listdir(campaign_dir):
        campaign_path = os.path.join(campaign_dir, campaign_name)
        campaign = pd.read_csv(campaign_path, keep_default_na=False)
        for _, row in tqdm(campaign.iterrows(), desc='analyzing variants'):
            enzyme = row['enzyme']
            substrate = row['substrate']
            var = row['var']
            var_folder_format = var.replace(':', '/')
            struc_pdb_path = os.path.join(structure_dir, f'{enzyme}-{substrate}', f'{var}_0.pdb')

            report_path = os.path.join(results_dir, f'{enzyme}-{substrate}_{var}')
            working_directory = '/disk2/lukas/plip'
            if not os.path.exists(report_path):
                os.makedirs(report_path)
            cli = ['python',
                '-m',
                'plip.plipcmd',
                '-f',
                struc_pdb_path,
                '--out',
                report_path,
                '--xml']
            cli_return = subprocess.run(cli, cwd=working_directory)
            #print(cli_return)


# Create a ZS fitness comparison df
create_comparison_df_flag = False
comparison_df_path = os.path.join(results_dir, f'comparison_df.csv')
if create_comparison_df_flag:
    comparison_df = pd.DataFrame(columns=['campaign_name', 'mutated_positions', 'active_site_plip', 'active_site_5A', 'active_site_10A', 'fitness', 'ZS1', 'ZS2', 'ZS3', 'ZS4', 'ZS5', 'ZS6', 'ZS7', 'ZS8', 'ZS9', 'ZS10'])
    # add all variants to the ZS fitness comparison df
    for campaign_path in dirpaths(campaign_dir): 
        campaign = pd.read_csv(campaign_path)
        for i in tqdm(range(0, campaign.shape[0]), desc=f'Creating comparison df, reading in campaign {os.path.basename(campaign_path)}...'):
            variant = campaign.iloc[i, :]
            # add active site to each variant
            mutated_pos = extract_mutated_residues(campaign, i)
            var = variant['var']
            xml_path = os.path.join(results_dir, f'{campaign_name}_{var}', 'report.xml')
            active_site_res = extract_active_site_from_xml(xml_path)
            pdb_path = os.path.join(structure_dir, campaign_name, f'{var}_0.pdb')
            active_site_5A = extract_active_site_by_radius(pdb_path, 5)
            active_site_10A = extract_active_site_by_radius(pdb_path, 10)
            new_row = pd.DataFrame({'campaign_name': campaign_name, 'mutated_positions':[campaign.iloc[i]['var']], 'active_site_plip':[active_site_res], 'active_site_5A':[active_site_5A], 'active_site_10A':[active_site_10A], 'fitness':[campaign.iloc[i]['fitness']], 'ZS1':[0.], 'ZS2':[0.], 'ZS3':[0.], 'ZS4':[0.], 'ZS5':[0.], 'ZS6':[0.], 'ZS7':[0.], 'ZS8':[0.], 'ZS9':[0.], 'ZS10':[0.], })
            comparison_df = pd.concat([comparison_df, new_row], ignore_index=True)
    # Save the comparison df
    comparison_df.to_csv(comparison_df_path, index=False)
    print(f'Saved comparison df to {comparison_df_path}')
else:
    comparison_df = pd.read_csv(comparison_df_path)


# iterate over all variants in comparison df 
for n, row in tqdm(comparison_df.iterrows(), desc='Calculating ZS for variants..'):
    mutated_pos = row['mutated_positions']
    mutated_pos_idx = re.findall(r'\d+', mutated_pos)
    new_res = [mutation[-1] for mutation in mutated_pos.split(':')]
    campaign_name = row['campaign_name']

    # Hydropathy of active site (kyte-doolittle) for plip active site
    hydropathy_score = 0
    active_site_res = row['active_site_plip']
    if isinstance(active_site_res, str):
        active_site_res = ast.literal_eval(active_site_res)
    for res in active_site_res.values():
        hydropathy_score += kyte_doolittle[res]
    hydropathy_score = hydropathy_score / len(active_site_res) 
    # normalizing the kyte-doolittle hydropathy score [0:1]
    kyte = (hydropathy_score - (min(kyte_doolittle.values()))) / (max(kyte_doolittle.values()) - (min(kyte_doolittle.values())))

    # Hydropathy of ligands (octanol-water partition) for plip active site
    struc_pdb_path = os.path.join(structure_dir, f'{campaign_name}', f'{mutated_pos}_0.pdb')
    mol = extract_ligands_to_mol(struc_pdb_path)
    logp = Descriptors.MolLogP(mol)
    # normalizing the octanol-water partition coefficient [0:1]
    min_logp = -3
    max_logp = 10
    logp = (logp - (min_logp)) / (max_logp - (min_logp))             


    # Hydropathy of active site (kyte-doolittle) for 5A active site
    hydropathy_score_5A = 0
    active_site_res = row['active_site_5A']
    if isinstance(active_site_res, str):    
        active_site_res = ast.literal_eval(active_site_res)
    for res in active_site_res.values():
        hydropathy_score_5A += kyte_doolittle[res]
    hydropathy_score_5A = hydropathy_score_5A / len(active_site_res) 
    # normalizing the kyte-doolittle hydropathy score [0:1]
    kyte_5A = (hydropathy_score_5A - (min(kyte_doolittle.values()))) / (max(kyte_doolittle.values()) - (min(kyte_doolittle.values())))


    # Hydropathy of active site (kyte-doolittle) for 5A active site
    hydropathy_score_10A = 0
    active_site_res = row['active_site_10A']
    if isinstance(active_site_res, str):
        active_site_res = ast.literal_eval(active_site_res)
    for res in active_site_res.values():
        hydropathy_score_10A += kyte_doolittle[res]
    hydropathy_score_10A = hydropathy_score_10A / len(active_site_res) 
    # normalizing the kyte-doolittle hydropathy score [0:1]
    kyte_10A = (hydropathy_score_10A - (min(kyte_doolittle.values()))) / (max(kyte_doolittle.values()) - (min(kyte_doolittle.values())))


    # plip report df
    xml_path = os.path.join(results_dir, f'{campaign_name}_{mutated_pos}', 'report.xml')
    tree = ET.parse(xml_path)
    root = tree.getroot()

    interaction_types = {
    'hydrophobic_interactions': ['dist'],
    'hydrogen_bonds': ['dist_h-a', 'dist_d-a', 'don_angle'],
    'pi_stacks': ['centdist', 'angle', 'offset'],
    'salt_bridges': ['dist'],  # Removed 'protispos' from fields
    'pi_cation_interactions': ['dist', 'offset'], 
    'water_bridges': ['dist_a-w', 'dist_d-w', 'water_angle'],
    'halogen_bonds': ['dist', 'don_angle', 'acc_angle'],
    'metal_complexes': ['dist', 'rms']
    }

    data = []
    for interaction_type, fields in interaction_types.items():
        parent = root.find(f".//interactions/{interaction_type}")
        if parent is not None:
            for interaction in parent:
                row = {'interaction_type': interaction_type}
                for field in fields:
                    element = interaction.find(field)
                    if element is not None and element.text:
                        text = element.text.strip()
                        if is_number(text):
                            row[field] = float(text)
                        elif text.lower() in ['true', 'false']:
                            row[field] = text.lower() == 'true'
                        else:
                            row[field] = text
                    else:
                        row[field] = np.nan
                data.append(row)
    plip_df = pd.DataFrame(data)
    #plip_df.to_csv('/disk2/lukas/tmp.csv')
    
    # Plip metrics 
    n_hydrophobic = plip_df[plip_df['interaction_type'] == 'hydrophobic_interactions'].shape[0]
    n_hydrogen_bonds = plip_df[plip_df['interaction_type'] == 'hydrogen_bonds'].shape[0]
    n_pi_stacks = plip_df[plip_df['interaction_type'] == 'pi_stacks'].shape[0]
    n_salt_bridges = plip_df[plip_df['interaction_type'] == 'salt_bridges'].shape[0]
    n_pi_cation_interactions = plip_df[plip_df['interaction_type'] == 'pi_cation_interactions'].shape[0]
    n_water_bridges = plip_df[plip_df['interaction_type'] == 'water_bridges'].shape[0]
    n_halogen_bonds = plip_df[plip_df['interaction_type'] == 'halogen_bonds'].shape[0]
    n_metal_complexes = plip_df[plip_df['interaction_type'] == 'metal_complexes'].shape[0]

    interaction_types = plip_df['interaction_type'].unique()
    averages = {}

    for interaction in interaction_types:
        subset = plip_df[plip_df['interaction_type'] == interaction]
        numeric_columns = subset.select_dtypes(include=[float, int]).columns
        averages[interaction] = subset[numeric_columns].apply(
            lambda col: col[col > 0].mean() if col.name != 'interaction_type' else None
        ).dropna().to_dict()


    # ZS1 only active site
    ZS1 = kyte
    ZS1_desc = 'Kyte'

    # ZS2 only active site
    ZS2 =  logp
    ZS2_desc = 'logp'

    # ZS3 sum active site and ligand (should behave same as delta)
    ZS3 = kyte + logp
    ZS3_desc = 'Kyte + logp'

    # ZS4 sum active site and ligand penalized by discrepancy of them
    ZS4 = (kyte + logp) * (1 - abs(kyte - logp))
    ZS4_desc = '(kyte + logp) * (1 - abs(kyte - logp))'

    # ZS5 only the delta between active site and ligand (no abs())
    ZS5 = kyte - logp
    ZS5_desc = 'Kyte - logp'    


    # ZS6
    ZS6 = pow((n_hydrophobic + n_hydrogen_bonds + n_pi_stacks + n_salt_bridges + n_pi_cation_interactions + n_water_bridges + n_pi_cation_interactions + n_water_bridges + n_halogen_bonds + n_metal_complexes), (kyte - logp))
    ZS6_desc = 'pow(n_iteraction, (kyte - logp))'

    # ZS7
    ZS7 = n_hydrophobic + n_hydrogen_bonds + n_pi_stacks + n_salt_bridges + n_pi_cation_interactions + n_water_bridges + n_pi_cation_interactions + n_water_bridges + n_halogen_bonds + n_metal_complexes
    ZS7_desc = 'n_interactions from plip'

    # ZS8 
    ZS8 = (
    (n_hydrophobic * averages.get('hydrophobic_interactions', {}).get('dist', 0)) +
    (n_hydrogen_bonds * ((averages.get('hydrogen_bonds', {}).get('dist_h-a', 0) + averages.get('hydrogen_bonds', {}).get('dist_d-a', 0)) / 2)) +
    (n_pi_stacks * averages.get('pi_stacks', {}).get('centdist', 0)) +
    (n_salt_bridges * averages.get('salt_bridges', {}).get('dist', 0)) +
    (n_pi_cation_interactions * averages.get('pi_cation_interactions', {}).get('dist', 0)) +
    (n_water_bridges * averages.get('water_bridges', {}).get('dist', 0)) +
    (n_halogen_bonds * averages.get('halogen_bonds', {}).get('dist', 0)) +
    (n_metal_complexes * averages.get('metal_complexes', {}).get('dist', 0))
    )
    ZS8_desc = 'Just plip interactions, more granular)'

    ZS9 = kyte_5A - logp 
    ZS9_desc = 'kyte_5A - logp'


    ZS10 = kyte_10A - logp
    ZS10_desc = 'kyte_10A - logp'
    
    # add the ZS to the comparison df
    comparison_df.at[n, 'ZS1'] = ZS1
    comparison_df.at[n, 'ZS2'] = ZS2
    comparison_df.at[n, 'ZS3'] = ZS3
    comparison_df.at[n, 'ZS4'] = ZS4
    comparison_df.at[n, 'ZS5'] = ZS5
    comparison_df.at[n, 'ZS6'] = ZS6
    comparison_df.at[n, 'ZS7'] = ZS7
    comparison_df.at[n, 'ZS8'] = ZS8
    comparison_df.at[n, 'ZS9'] = ZS9
    comparison_df.at[n, 'ZS10'] = ZS10


gt = comparison_df['fitness'].to_numpy()
ZS1 = comparison_df['ZS1'].to_numpy()
ZS2 = comparison_df['ZS2'].to_numpy()
ZS3 = comparison_df['ZS3'].to_numpy()
ZS4 = comparison_df['ZS4'].to_numpy()
ZS5 = comparison_df['ZS5'].to_numpy()
ZS6 = comparison_df['ZS6'].to_numpy()
ZS7 = comparison_df['ZS7'].to_numpy()
ZS8 = comparison_df['ZS8'].to_numpy()
ZS9 = comparison_df['ZS9'].to_numpy()
ZS10 = comparison_df['ZS10'].to_numpy()

rho1, _ = stats.spearmanr(gt, ZS1)
rho2, _ = stats.spearmanr(gt, ZS2)
rho3, _ = stats.spearmanr(gt, ZS3)
rho4, _ = stats.spearmanr(gt, ZS4)
rho5, _ = stats.spearmanr(gt, ZS5)
rho6, _ = stats.spearmanr(gt, ZS6)
rho7, _ = stats.spearmanr(gt, ZS7)
rho8, _ = stats.spearmanr(gt, ZS8)
rho9, _ = stats.spearmanr(gt, ZS9)
rho10, _ = stats.spearmanr(gt, ZS10)


print(f'{ZS1_desc}: {rho1} \n {ZS2_desc}: {rho2} \n {ZS3_desc}: {rho3} \n {ZS4_desc}: {rho4} \n {ZS5_desc}: {rho5} \n {ZS6_desc}: {rho6} \n {ZS7_desc}: {rho7} \n {ZS8_desc}: {rho8} \n {ZS9_desc}: {rho9} \n {ZS10_desc}: {rho10} ')



ipdb.set_trace()
