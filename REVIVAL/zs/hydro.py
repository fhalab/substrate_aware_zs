# File to automatically run hydrophobicity calculations as a ZS
# Author: Lukas Radtke {lradtke@caltech.edu}
# Date: 16/12/24

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
from scipy import stats
import ipdb
from rdkit import Chem
from openbabel import pybel

from REVIVAL.preprocess import ZSData

class HydroData(ZSData):
    def __init__(
            self,
            input_csv = str,
            scale_fit: str = 'parent',
            combo_col_name: str = 'AAs',
            var_col_name: str = 'var',
            fit_col_name: str = 'fitness',
            enzyme_col_name: str = 'enzyme',
            seq_col_name: str = 'seq',
            substrate_col_name: str = 'substrate',
            cofactor_col_name: str = 'cofactor',
            substrate_smiles_col_name: str = 'substrate-smiles',
            cofactor_smiles_col_name: str = 'cofactor-smiles',
            
            structure_dir: str = 'data/structure',
            plip_dir: str = '../plip',
            results_dir: str = 'results',   
            active_site_radius: int = 10, 
            variant_structures_available: bool = True,
            variant_structure_dir: str = '/disk2/lukas/EnzymeOracle/data/struc_test/struc/chai_run_0'
        ):

        super().__init__(
            input_csv = input_csv,
            scale_fit = scale_fit,
            combo_col_name = combo_col_name,
            var_col_name = var_col_name,
            fit_col_name = fit_col_name,
            structure_dir = structure_dir,
        )

        self.plip_dir = plip_dir
        self.enzyme_col_name = enzyme_col_name
        self.substrate_col_name = substrate_col_name
        self.seq_col_name = seq_col_name
        self.cofactor_col_name = cofactor_col_name
        self.substrate_smiles_col_name = substrate_smiles_col_name
        self.cofactor_smiles_col_name = cofactor_smiles_col_name
        self.structure_dir = structure_dir
        self.results_dir = results_dir 
        self.active_site_radius = active_site_radius
        self.variant_structures_available = variant_structures_available    
        self.variant_structure_dir = variant_structure_dir

        self.hopp_woods_scale = {
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

        self.oneletter2threeletter =  {
                                    'A': 'ALA',
                                    'R': 'ARG', 
                                    'N': 'ASN',  
                                    'D': 'ASP',  
                                    'C': 'CYS', 
                                    'E': 'GLU', 
                                    'Q': 'GLN',  
                                    'G': 'GLY',
                                    'H': 'HIS',  
                                    'I': 'ILE',  
                                    'L': 'LEU',
                                    'K': 'LYS',
                                    'M': 'MET',  
                                    'F': 'PHE',  
                                    'P': 'PRO',
                                    'S': 'SER',
                                    'T': 'THR',  
                                    'W': 'TRP',  
                                    'Y': 'TYR', 
                                    'V': 'VAL',
                                }

        self.inference()


    def extract_mutated_residues(self, csv):
        """
        Extracts the mutations position and prepares them in the expected format for hydro (e.g. 'A165 A183 A301').
        Assumes the same mutation positions within the entire campaign.
        """
        position_list = []
        for _, row in csv.iterrows():
            var = row[self._var_col_name]
            mutated_pos = re.findall(r"\d+", var)
            all_mutations = ['A'+mut for mut in mutated_pos]
            for res in all_mutations:
                if res not in position_list:
                    position_list.append(res)
        mutated_pos = ' '.join(position_list)
        return mutated_pos
    

    def extract_ca_coordinates(self, pdb_file):
        """
        Extracts the coordinates of C_alpha atoms from a pdb file and returns them as a dict.
        The dict is of format: key=residue_position_in_sequence: value=cartesian coordinates of C_alpha atom
        """
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


    def extract_hetatm_coordinates(self, pdb_file):
        """
        Extracts all heteroatoms from a pdb file. 
        Returns a nested list of all cartesian coordinates.
        """
        coordinates = []
        with open(pdb_file, 'r') as file:
            for line in file:
                if line.startswith('HETATM') or 'LIG' in line:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        coordinates.append([x, y, z])
        return coordinates


    def extract_active_site_by_radius_parent(self, pdb_path, radius_angstrom, mut_pos_idx, variant):
        """
        Extracts the active site residues from the parent PDB file using a radius from the ligand in cartesian space.
        Residues with their C-alpha atom within the radius will be added to the return dictionar with format: 
        dict {key=str(index): value=str(threelettercode), }
        """
        hetatm_coords = self.extract_hetatm_coordinates(pdb_path)
        ca_coords_idx = self.extract_ca_coordinates(pdb_path)
        ca_coords = np.array(list(ca_coords_idx.values()))
        ca_indices = np.array(list(ca_coords_idx.keys()))
        active_site_ca_idx = []
        for hetatm_coord in hetatm_coords:
            hetatm_coords = np.array(hetatm_coord)
            indices_within_radius = np.where(np.linalg.norm(ca_coords - hetatm_coord, axis=1) <= radius_angstrom)[0]
            active_site_ca_idx.extend(ca_indices[indices_within_radius])
        active_site_ca_idx = set(active_site_ca_idx)
        active_site = {}
        for active_site_idx in active_site_ca_idx:
            active_site[str(active_site_idx)] = variant[self.seq_col_name][active_site_idx]
        for i, mut_pos in enumerate(mut_pos_idx):
            active_site[mut_pos] = variant[self._combo_col_name][i]     
        return active_site


    def extract_ligands_to_mol(self, pdb):
        """
        Extracts the ligand from a PDB file and returns it as a mol object.
        """
        hetatm_lines = []
        with open(pdb, 'r') as file:
            for line in file:
                if line.startswith('ATOM') and 'LIG' in line or line.startswith('HETATM'):
                    hetatm_lines.append(line)
        with open('temp_hetatm.pdb', 'w') as temp_file:
            temp_file.writelines(hetatm_lines)
        mol = next(pybel.readfile('pdb', 'temp_hetatm.pdb'))
        if mol == None:
            print('ERROR: could not parse molecule from PDB file.')
            ipdb.set_trace()
        else:
            try: 
                Chem.SanitizeMol(mol)
            except Exception as e:
                print(f'Sanitation failed: {e}')
        return mol
    

    def extract_active_site_from_xml(self, xml_path):
        """
        Reads a .xml file generated by PLIP and extracts the active site residues.
        Return active site dictionary of format: 
        """
        with open(xml_path, 'r') as xml:
            report = xmltodict.parse(xml.read())
        bs_residues = report['report']['bindingsite']['bs_residues']['bs_residue']
        active_site_res_parent = {} #position: identity
        for bs_residue in bs_residues:
            if bs_residue['@contact'] == 'False':
                active_site_res_parent[bs_residue['#text'][:-1]] = bs_residue['@aa']
        return active_site_res_parent

    
    def run_plip(self, csv):
        """
        Runs PLIP for each variant of a campaign csv and saves the plip report as .xml file.
        If there is a PLIP report already present, the function omits creating a new one.
        If a pdb file for the respective variant is not available, the parent structure will be used as a fallback.
        """
        for _, row in tqdm(csv.iterrows(), desc='Analyzing variants with PLIP...', ncols=1000):
            enzyme = row[self.enzyme_col_name]
            substrate = row[self.substrate_col_name]
            var = row[self._var_col_name]
            struc_pdb_path = os.path.join(self.variant_structure_dir, f'{enzyme}-{substrate}', f'{var}_0.pdb')

            report_path = os.path.join(self.results_dir, 'plip_reports', f'{enzyme}-{substrate}_{var}')
            
            # check if a plip report has already been saved for this variant.
            check_report_path = os.path.join(report_path, 'report.xml')
            #if check_report_path == '/disk2/lukas/EnzymeOracle/data/results/plip_reports/PfTrpB-5cyano_I165K:I183G:Y301K/report.xml':
            #    ipdb.set_trace()
            if os.path.exists(check_report_path) == True:
                pass
            elif os.path.exists(struc_pdb_path) == False:
                print(f'ERROR: Did not find pdb file for variant {os.path.basename(struc_pdb_path)}, using parent structure for PLIP.')
                struc_pdb_path = os.path.join(self.structure_dir, f'{enzyme}-{substrate}.pdb')
                working_directory = self.plip_dir
                os.makedirs(report_path, exist_ok=True)
                cli = ['python',
                    '-m',
                    'plip.plipcmd',
                    '-f',
                    struc_pdb_path,
                    '--out',
                    report_path,
                    '--xml']
                subprocess.run(cli, cwd=working_directory)
            else:
                working_directory = self.plip_dir
                os.makedirs(report_path, exist_ok=True)
                cli = ['python',
                    '-m',
                    'plip.plipcmd',
                    '-f',
                    struc_pdb_path,
                    '--out',
                    report_path,
                    '--xml']
                subprocess.run(cli, cwd=working_directory)
                if os.path.exists(check_report_path) == False and check_report_path == '/disk2/lukas/EnzymeOracle/data/results/plip_reports/PfTrpB-5cyano_I165K:I183G:Y301K/report.xml':
                    ipdb.set_trace()
    

    def calc_hydrophobitcity(self, active_site, scale):
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
                hydrophobicity += scale[self.oneletter2threeletter[res]]
        hydrophobicity = hydrophobicity / len(active_site) 
        return hydrophobicity


    def inference(self):
        # read the campaign data
        csv = pd.read_csv(self._input_csv, keep_default_na=False)
        csv = csv[~csv.apply(lambda row: row.astype(str).str.contains(r"\*").any(), axis=1)].reset_index(drop=True)  # filters for stop codons, marked as '*'
        
        # load parend pdb for campaign
        enzyme_name = csv[self.enzyme_col_name][1] # assumes its the same enzyme for every variant
        substrate_name = csv[self.substrate_col_name][1] # assumes its the same enzyme for every variant
        pdb = os.path.join(self.structure_dir, f'{enzyme_name}-{substrate_name}.pdb')

        # prepares mutated position for the radius-based active site definition
        campaign_name = os.path.basename(self._input_csv[:-4])
        mutated_positions = self.extract_mutated_residues(csv)
        mut_pos_idx = [pos[1:] for pos in mutated_positions.split()]

        if self.variant_structures_available == True:
            self.run_plip(csv)
        

        # initializing score_df with active sites 
        rows = []
        for _, variant in tqdm(csv.iterrows(), desc=f'Creating comparison df, reading in campaign {os.path.basename(self._input_csv)}...'):
            
            # generating radius-based active site 
            active_site_radius = self.extract_active_site_by_radius_parent(pdb, self.active_site_radius, mut_pos_idx, variant)
            
            # generating plip-baed active site 
            if self.variant_structures_available == True:
                xml_path = os.path.join(self.results_dir, 'plip_reports', f'{variant[self.enzyme_col_name]}-{variant[self.substrate_col_name]}_{variant[self._var_col_name]}', 'report.xml')
                active_site_plip = self.extract_active_site_from_xml(xml_path)
                for i, mut_pos in enumerate(mut_pos_idx):
                    active_site_plip[mut_pos] = variant[self._combo_col_name][i]  
            else:
                active_site_plip = {}


            gt_fitness = variant[self._fit_col_name]
            new_row = {'campaign_name': campaign_name, 'mutated_positions':[mutated_positions], 'active_site_radius':[active_site_radius], 'active_site_plip':[active_site_plip], 'gt_fitness':gt_fitness, 'ZS1':[0.], 'ZS2':[0.]}
            rows.append(new_row)

        score_df = pd.DataFrame(data=rows)

        mol = self.extract_ligands_to_mol(pdb)


        # add the two ZS scores to the score_df
        for n, row in tqdm(score_df.iterrows(), desc=f'Calulating ZS for variants...'):
            
            # calculate hydrophobicity with hopp scale and 10A active site definition
            hydro_hopp_radius = self.calc_hydrophobitcity(row['active_site_radius'], self.hopp_woods_scale)
            hydro_hopp_radius = (hydro_hopp_radius - (min(self.hopp_woods_scale.values()))) / (max(self.hopp_woods_scale.values()) - (min(self.hopp_woods_scale.values())))
            
            # calculate hydrophobicity with hopp scale and PLIP active site definition
            if self.variant_structures_available == True:
                hydro_hopp_PLIP = self.calc_hydrophobitcity(row['active_site_plip'], self.hopp_woods_scale)
            else:
                pass

            # calculate hydrophobicity (tpsa) of ligand
            tpsa = mol.calcdesc()["TPSA"]

            # calculate the ZS for radius-based active site definition
            ZS1 = tpsa - hydro_hopp_radius
            score_df.at[n, 'ZS1'] = ZS1

            # calculate the ZS for plip-based active site definition
            if self.variant_structures_available == True:
                ZS2 = tpsa - hydro_hopp_PLIP
            else:
                ZS2 = 0
            score_df.at[n, 'ZS2'] = ZS2

        """
        gt = score_df['gt_fitness'].to_numpy()
        ZS1 = score_df['ZS1'].to_numpy()
        ZS2 = score_df['ZS2'].to_numpy()

        rho1, _ = stats.spearmanr(gt, ZS1)
        rho2, _ = stats.spearmanr(gt, ZS2)
        """

        # save the score df
        score_df_path = os.path.join(self.results_dir, f'{os.path.basename(self._input_csv)[:-4]}_score_df.csv')
        score_df.to_csv(score_df_path)
        print(f'ZS score .csv file saved to {score_df_path}.')


def run_hydro(pattern: str | list = None, kwargs: dict = {}):    
    
    if isinstance(pattern, str):
        lib_list = glob(pattern)
    else:
        lib_list = deepcopy(pattern)
    
    for lib in lib_list:
        HydroData(input_csv=lib, **kwargs)


run_hydro(pattern='/disk2/lukas/EnzymeOracle/data/multi_substrate/meta/scale2parent/*.csv', 
            kwargs={'active_site_radius': 15,
                'variant_structures_available': True,
                'structure_dir': '../EnzymeOracle/data/multi_substrate/structure',
                'results_dir': '/disk2/lukas/EnzymeOracle/data/results',
                }
            )



