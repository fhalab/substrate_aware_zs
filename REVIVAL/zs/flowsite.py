# File to automatically run Flowsite as a ZS
# Author: Lukas Radtke {lradtke@caltech.edu}
# Date: 06/12/24

import numpy as np
import pandas as pd
import os 
import subprocess
import re
import ipdb
from glob import glob
from copy import deepcopy
from tqdm import tqdm
from scipy import stats
from Bio.PDB import PDBParser

from REVIVAL.preprocess import ZSData
from REVIVAL.util import checkNgen_folder

class FlowsiteData(ZSData):
    def __init__(
            self,
            input_csv = str,
            scale_fit: str = 'parent',
            structure_dir: str = 'data/structure',
            flowsite_dir: str = '../FlowSite',
            results_dir: str = './outputs',
            combo_col_name: str = 'AAs',
            var_col_name: str = 'var',
            seq_col_name: str = 'seq',
            fit_col_name: str = 'fitness',
            substrate_col_name: str = 'substrate',
            cofactor_col_name: str = 'cofactor',
            substrate_smiles_col_name: str = 'substrate-smiles',
            cofactor_smiles_col_name: str = 'cofactor-smiles',
            num_inference: int = 1,
            batch_size: int = 10,
            fake_constant_dur: int = 100000,
            fake_decay_dur: int = 10,
            fake_ratio_start: float = 0.2,
            fake_ratio_end: float = 0.2,
            residue_loss_weight: float = 0.2,
            flow_matching_sigma: float = 0.5,
            prior_scale: int = 1,
            num_angle_pred: int = 5,
            pocket_residue_cutoff: int = 12,

    ):
        
        super().__init__(
            input_csv = input_csv,
            scale_fit = scale_fit,
            combo_col_name = combo_col_name,
            var_col_name = var_col_name,
            seq_col_name = seq_col_name,
            fit_col_name = fit_col_name,
            structure_dir = structure_dir,
        )
          
        self.flowsite_dir = flowsite_dir
        self.inference_path = os.path.join(self.flowsite_dir, 'inference.py')
        self.substrate_col_name = substrate_col_name
        self.cofactor_col_name = cofactor_col_name
        self.substrate_smiles_col_name = substrate_smiles_col_name
        self.cofactor_smiles_col_name = cofactor_smiles_col_name
        self.structure_dir = structure_dir
        self.flowsite_dir = flowsite_dir
        self.results_dir = results_dir     
        self.num_inference = num_inference
        self.batch_size = batch_size
        
        self.model_dict = {1: 'pocket_gen/lf5t55w4/checkpoints/best.ckpt',
                            2: 'pocket_gen/b1ribx1a/checkpoints/best.ckpt'}
        self.model = 2     # recommended to use model 2 when using --pocked_def_residues
        if self.model == 1:
            self.ns = 48
            self.nv = 12
        elif self.model == 2:
            self.ns = 32
            self.nv = 8

        self.fake_constant_dur = fake_constant_dur
        self.fake_decay_dur = fake_decay_dur
        self.fake_ratio_start = fake_ratio_start
        self.fake_ratio_end = fake_ratio_end
        self.residue_loss_weight = residue_loss_weight
        self.flow_matching_sigma = flow_matching_sigma
        self.prior_scale = prior_scale
        self.num_angle_pred = num_angle_pred
        self.pocket_residue_cutoff = pocket_residue_cutoff

        print('Running flowsite...')
        self.inference2()

        print(f"Creating score .csv for {self.lib_name}")
        self.score_mut_file()


    def extract_design_residues(self, row):
        var = row[self._var_col_name]
        mutated_pos = re.findall(r"\d+", var)
        all_mutations = ['A'+mut for mut in mutated_pos]
        design_residues = ' '.join(all_mutations)    
        return design_residues
    

    def extract_smiles(self, row):
        substrate_smiles = row[self.substrate_smiles_col_name]
        if self.cofactor_col_name in row.keys():
            cofactor_smiles = row[self.cofactor_smiles_col_name]
            return substrate_smiles+'.'+cofactor_smiles
        else:
            return substrate_smiles
    
    def extract_muation_positions(self, csv):
        position_list = []
        for _, row in csv.iterrows():
            var = row[self._var_col_name]
            mutated_pos = re.findall(r"\d+", var)
            all_mutations = ['A'+mut for mut in mutated_pos]
            for res in all_mutations:
                if res not in position_list:
                    position_list.append(res)
        return position_list


    def inference2(self) -> None:
        """
        Prepares the inputs for the model and runs FlowSite
        """
        csv = pd.read_csv(self._input_csv, keep_default_na=False)
        csv = csv[~csv.apply(lambda row: row.astype(str).str.contains(r"\*").any(), axis=1)].reset_index(drop=True)  # filters for stop codons, marked as '*'
        
        # find the name of the enzyme and substrate in the campaign
        enzyme_name = csv['enzyme'][1] # assumes its the same enzyme for every variant
        substrate_name = csv['substrate'][1] # assumes its the same enzyme for every variant

        # find enzyme structure file based on enzyme name and ligand name
        pdb = os.path.join(self.structure_dir, f'{enzyme_name}-{substrate_name}.pdb')
        assert os.path.exists(pdb), f'No pdb file found at path {pdb}, structure name should follow the format: <enzyme_name>-<substrate_name>.pdb'
        
        # create a smiles for substrate and possibly cofactor
        # find ligand.sdf file, possibly with cofactor included  
        substrate_smiles = csv['substrate-smiles'][1]
        if self.cofactor_col_name in csv.columns:
            cofactor_smiles = csv['cofactor-smiles'][1]
            smiles = substrate_smiles+'.'+cofactor_smiles
            cofactor_name = csv['cofactor'][1]
            ligand_sdf = os.path.join(self.structure_dir, f'{substrate_name}_{cofactor_name}.sdf')
        else:
            smiles = substrate_smiles
            ligand_sdf = os.path.join(self.structure_dir, f'{substrate_name}.sdf')

        # extract all mutated positions from the campaign to inform for which positions to generate for
        mutated_pos_list = self.extract_muation_positions(csv)
        mutated_pos_list = ', '.join(mutated_pos_list)

        # extract chain length
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('PDB_structure', pdb)
        residues = [res for res in structure.get_residues() if res.id[0] == ' ']
        num_residues = len(residues)

        design_all_str = f'A0-{num_residues}'

        # load model
        model_path = os.path.join(self.flowsite_dir, self.model_dict[self.model] )

        # define output file name
        output_dir = os.path.join(self.results_dir, f'{enzyme_name}-{substrate_name}_inference')
        os.mkdir(output_dir, exist_ok=True)

        env = os.environ.copy()
        env['CUDA_VISIBLE_DEVICES'] = '0'

        cli_cmd = [ 'python', self.inference_path,
                    '--num_inference', str(self.num_inference),
                    '--batch_size', str(self.batch_size),
                    '--out_dir', self.dir,
                    '--protein', pdb,
                    '--pocket_def_residues', mutated_pos_list,
                    '--smiles', smiles,
                    '--design_residues', mutated_pos_list,
                    '--checkpoint', model_path,
                    '--run_test', 
                    '--run_name', 'inference1',
                    '--layer_norm', 
                    '--fake_constant_dur', str(self.fake_constant_dur),
                    '--fake_decay_dur', str(self.fake_decay_dur),
                    '--fake_ratio_start', str(self.fake_ratio_start),
                    '--fake_ratio_end', str(self.fake_ratio_end),
                    '--residue_loss_weight', str(self.residue_loss_weight),
                    '--use_tfn',
                    '--use_inv',
                    '--time_condition_tfn',
                    '--correct_time_condition',
                    '--time_condition_inv',
                    '--time_condition_repeat',
                    '--flow_matching',
                    '--flow_matching_sigma', str(self.flow_matching_sigma),
                    '--prior_scale', str(self.prior_scale),
                    '--num_angle_pred', str(self.num_angle_pred),
                    '--ns', str(self.ns),
                    '--nv', str(self.nv),
                    '--batch_norm',
                    '--self_condition_x',
                    '--self_condition_inv',
                    '--no_tfn_self_condition_inv',
                    '--self_fancy_init',
                    '--pocket_residue_cutoff', str(self.pocket_residue_cutoff)
        ]

        cli_return = subprocess.run(cli_cmd, capture_output=True, text=True, env=env)
        print(cli_return)
        ipdb.set_trace()


    def score_mut_file(self) -> None:
        """
        A function for extracting the ZS scores from Flowsite output files and writing them to a csv.
        """
        score_df = pd.DataFrame(columns=["variant", "gt_fitness", "zs_score"])
        csv = pd.read(self._input_csv, keep_default_na=False)
        csv = csv[~csv.apply(lambda row: row.astype(str).str.contains(r"\*").any(), axis=1)].reset_index(drop=True)  # filters for stop codons, marked as '*'

        enzyme_name = csv['enzyme'][1] # assumes its the same enzyme for every variant
        substrate_name = csv['substrate'][1] # assumes its the same enzyme for every variant
        output_dir = os.path.join(self.results_dir, f'{enzyme_name}-{substrate_name}_inference') 

        design_logits = np.load('_') 
        ############## RESUME HERE



    def inference(self) -> None:
        """
        Prepares the inputs for the model and runs FlowSite
        """
        csv = pd.read_csv(self._input_csv, keep_default_na=False)
        csv = csv[~csv.apply(lambda row: row.astype(str).str.contains(r"\*").any(), axis=1)].reset_index(drop=True)  # filters for stop codons, marked as '*'
        
        # create inference  with cols smiles, protein, sign_residues, pocket_def_residues
        inference_df = pd.DataFrame(columns=['ligand', 'smiles', 'protein', 'design_residues', 'pocket_def_ligand', 'pocket_def_residues', 'pocket_def_center' ])
        for idx, row in tqdm(csv.iterrows(), desc=f'Creating inference df for campaign {os.path.basename(self._input_csv)}...'):
            smiles = self.extract_smiles(row)
            structure_path = os.path.join(self.structure_dir, os.path.basename(self._input_csv).replace('.csv', '.pdb'))
            design_residues = self.extract_design_residues(row)
            inference_instance = pd.DataFrame(data = {'ligand': pd.NA, 'smiles': str(smiles), 'protein':  str(structure_path), 'design_residues': str(design_residues), 'pocket_def_ligand':pd.NA, 'pocket_def_center': pd.NA,'pocket_def_residues': str(design_residues), }, index=[0])
            if design_residues != '':
                inference_df = pd.concat([inference_df, inference_instance], ignore_index=True)
            else:
                print(f'No active sie defined in campaing {os.path.basename(self._input_csv)} and row {idx}. Skipping this entry...')

        inference_df_path = os.path.join(self.results_dir, f'inference_df_{os.path.basename(self._input_csv)}')
        inference_df.to_csv(inference_df_path, index=False, na_rep='')

        model_path = os.path.join(self.flowsite_dir, self.model_dict[self.model])

        env = os.environ.copy()
        env['CUDA_VISIBLE_DEVICES'] = '0'

        cli_cmd = [ 'python', self.inference_path,
                    '--num_inference', str(self.num_inference),
                    '--batch_size', str(self.batch_size),
                    '--out_dir', self.results_dir,
                    '--csv_file', inference_df_path,
                    '--checkpoint', model_path,
                    '--run_test', 
                    '--run_name', 'inference1',
                    '--layer_norm', 
                    '--fake_constant_dur', str(self.fake_constant_dur),
                    '--fake_decay_dur', str(self.fake_decay_dur),
                    '--fake_ratio_start', str(self.fake_ratio_start),
                    '--fake_ratio_end', str(self.fake_ratio_end),
                    '--residue_loss_weight', str(self.residue_loss_weight),
                    '--use_tfn',
                    '--use_inv',
                    '--time_condition_tfn',
                    '--correct_time_condition',
                    '--time_condition_inv',
                    '--time_condition_repeat',
                    '--flow_matching',
                    '--flow_matching_sigma', str(self.flow_matching_sigma),
                    '--prior_scale', str(self.prior_scale),
                    '--num_angle_pred', str(self.num_angle_pred),
                    '--ns', str(self.ns),
                    '--nv', str(self.nv),
                    '--batch_norm',
                    '--self_condition_x',
                    '--self_condition_inv',
                    '--no_tfn_self_condition_inv',
                    '--self_fancy_init',
                    '--pocket_residue_cutoff', str(self.pocket_residue_cutoff)
        ]

        print('Running flowsite...')
        cli_return = subprocess.run(cli_cmd, capture_output=True, text=True, env=env)
        print(cli_return)
        ipdb.set_trace()
    
    @property
    def pdb_path(self) -> str:
        """
        Filepath to the corresponding .pdb file for a given .csv of same name.
        """
        return os.path.join(
            self._structure_dir,
            os.path.basename(self._input_csv).replace(".csv", ".pdb"),
        )


def run_flowsite(pattern: str | list = None, kwargs: dict = {}):

    if isinstance(pattern, str):
        lib_list = glob(pattern)
    else:
        lib_list = deepcopy(pattern)

    for lib in lib_list:
        print(f'Running Flowsite for {lib}...')
        FlowsiteData(input_csv=lib, **kwargs) 
        
run_flowsite(pattern='/disk2/lukas/EnzymeOracle/data/multi_substrate/meta/scale2parent/*.csv', 
            kwargs={
                'flowsite_dir': '/disk2/lukas/FlowSite',
                'structure_dir': '../EnzymeOracle/data/multi_substrate/structure',
                'results_dir': '/disk2/lukas/EnzymeOracle/data/results'
                }
            )


"""
WORKS:
CUDA_VISIBLE_DEVICES="" python /disk2/lukas/FlowSite/inference.py --num_inference 1 --batch_size 10 --out_dir /disk2/lukas/EnzymeOracle/data/results --csv_file /disk2/lukas/EnzymeOracle/data/results/inference_df_PfTrpB-5iodo.csv --checkpoint /disk2/lukas/FlowSite/pocket_gen/b1ribx1a/checkpoints/best.ckpt --run_test --run_name inference1 --layer_norm --fake_constant_dur 100000 --fake_decay_dur 10 --fake_ratio_start 0.2 --fake_ratio_end 0.2 --residue_loss_weight 0.2 --use_tfn --use_inv --time_condition_tfn --correct_time_condition --time_condition_inv --time_condition_repeat --flow_matching --flow_matching_sigma 0.5 --prior_scale 1 --num_angle_pred 5 --ns 32 --nv 8 --batch_norm --self_condition_x 
"""