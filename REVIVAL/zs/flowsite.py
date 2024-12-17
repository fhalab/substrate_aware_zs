# File to automatically run Flowsite as a ZS
# Author: Lukas Radtke {lradtke@caltech.edu}
# Date: 06/12/24

import numpy as np
import pandas as pd
import os 
import subprocess
import re
import torch
from glob import glob
from copy import deepcopy
from tqdm import tqdm

from REVIVAL.preprocess import ZSData

class FlowsiteData(ZSData):
    def __init__(
            self,
            input_csv = str,
            scale_fit: str = 'parent',
            combo_col_name: str = 'AAs',
            var_col_name: str = 'var',
            fit_col_name: str = 'fitness',
            enzyme_col_name: str = 'enzyme', 
            substrate_col_name: str = 'substrate',
            cofactor_col_name: str = 'cofactor',
            substrate_smiles_col_name: str = 'substrate-smiles',
            cofactor_smiles_col_name: str = 'cofactor-smiles',
            
            structure_dir: str = 'data/structure',
            flowsite_dir: str = '../FlowSite',
            results_dir: str = 'results',
            
            num_inference: int = 10,
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
            fit_col_name = fit_col_name,
            structure_dir = structure_dir,
        )
          
        self.flowsite_dir = flowsite_dir
        self.inference_path = os.path.join(self.flowsite_dir, 'inference.py')
        self.enzyme_col_name = enzyme_col_name,
        self.substrate_col_name = substrate_col_name
        self.cofactor_col_name = cofactor_col_name
        self.substrate_smiles_col_name = substrate_smiles_col_name
        self.cofactor_smiles_col_name = cofactor_smiles_col_name
        self.structure_dir = structure_dir
        self.results_dir = results_dir     
        self.num_inference = num_inference
        self.batch_size = batch_size
        
        # model and hyperparameter choices according to https://github.com/HannesStark/FlowSite
        self.model_dict = {1: 'pocket_gen/lf5t55w4/checkpoints/best.ckpt',
                            2: 'pocket_gen/b1ribx1a/checkpoints/best.ckpt'}
        self.model = 2     
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

        self.oneletter2index = {'A': 0, 
                                'R': 1, 
                                'N': 2,
                                'D': 3, 
                                'C': 4, 
                                'Q': 5, 
                                'E': 6, 
                                'G': 7, 
                                'H': 8, 
                                'I': 9, 
                                'L': 10, 
                                'K': 11, 
                                'M': 12, 
                                'F': 13, 
                                'P': 14, 
                                'S': 15, 
                                'T': 16, 
                                'W': 17, 
                                'Y': 18, 
                                'V': 19}

        print(f'Running flowsite for {os.path.basename(self.lib_name)}...')
        self.inference()

        print(f"Creating score .csv for {os.path.basename(self.lib_name)}...")
        self.score_mut_file()


    def extract_smiles(self, csv):
        """
        Extracts substrate smiles and adds cofactor smiles if available.
        Assumes the same substrate and cofactor within the entire campaign.
        """
        substrate_smiles = csv['substrate-smiles'][1]
        if self.cofactor_col_name in csv.columns:
            cofactor_smiles = csv['cofactor-smiles'][1]
            return substrate_smiles+'.'+cofactor_smiles
        else:
            return substrate_smiles


    def extract_mutation_positions(self, csv):
        """
        Extracts the mutations position and prepares them in the expected format of flowsite (e.g. 'A165, A183, A301').
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
        mutated_pos = ', '.join(position_list)
        return mutated_pos


    def inference(self) -> None:
        """
        Prepares the inputs for the model and runs FlowSite
        """
        # read the campaign data
        csv = pd.read_csv(self._input_csv, keep_default_na=False)
        csv = csv[~csv.apply(lambda row: row.astype(str).str.contains(r"\*").any(), axis=1)].reset_index(drop=True)  # filters for stop codons, marked as '*'
        
        # find the name of the enzyme and substrate in the campaign
        enzyme_name = csv[self.enzyme_col_name][1] # assumes its the same enzyme for every variant
        substrate_name = csv[self.substrate_col_name][1] # assumes its the same enzyme for every variant

        # find enzyme structure file based on enzyme name and ligand name
        pdb = os.path.join(self.structure_dir, f'{enzyme_name}-{substrate_name}.pdb')
        assert os.path.exists(pdb), f'No pdb file found at path {pdb}, structure name should follow the format: <enzyme_name>-<substrate_name>.pdb'

        # extract all mutated positions from the campaign to inform for which positions to generate likelyhoods for
        mutated_pos = self.extract_mutation_positions(csv)

        # create a smiles for substrate and possibly cofactor   
        smiles = self.extract_smiles(csv)

        # load model
        model_path = os.path.join(self.flowsite_dir, self.model_dict[self.model])

        # define output file name
        output_dir = os.path.join(self.results_dir, f'{enzyme_name}-{substrate_name}_inference')
        os.makedirs(output_dir, exist_ok=True)

        # set up run device
        cuda_available = torch.cuda.is_available()
        env = os.environ.copy()
        if cuda_available:
            env['CUDA_VISIBLE_DEVICES'] = '0' 
        else:
            env.pop('CUDA_VISIBLE_DEVICES', None)

        # define console command and execute subprocess
        cli_cmd = [ 'python', self.inference_path,
                    '--num_inference', str(self.num_inference),
                    '--batch_size', str(self.batch_size),
                    '--out_dir', output_dir,
                    '--protein', pdb,
                    '--pocket_def_residues', mutated_pos,
                    '--smiles', smiles,
                    '--design_residues', mutated_pos,
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
        subprocess.run(cli_cmd, capture_output=True, text=True, env=env)


    def score_mut_file(self) -> None:
        """
        A function for extracting the ZS scores from flowsite output files and writing them to a csv.
        """

        score_df = pd.DataFrame(columns=["variant", "gt_fitness", "zs_score"])
        csv = pd.read_csv(self._input_csv, keep_default_na=False)
        csv = csv[~csv.apply(lambda row: row.astype(str).str.contains(r"\*").any(), axis=1)].reset_index(drop=True)  # filters for stop codons, marked as '*'

        enzyme_name = csv['enzyme'][1] # assumes its the same enzyme for every variant
        substrate_name = csv['substrate'][1] # assumes its the same enzyme for every variant
        output_dir = os.path.join(self.results_dir, f'{enzyme_name}-{substrate_name}_inference', 'complexid0', 'designed_logits.npy') 

        # load flowsite output
        ZS = np.load(output_dir)

        # extract ZS scores from flowsite output
        rows = []
        for _, variant in tqdm(csv.iterrows(), desc=f'Extracing ZS scores for {self._input_csv}'):            
            log_likelyhoods_avg = []
            for i, res in enumerate(variant[self._combo_col_name]):
                log_score = sum([ZS[inference][i][self.oneletter2index[res]] for inference in range(self.num_inference)])
                log_likelyhoods_avg.append(log_score)
            log_zs_score = sum(log_likelyhoods_avg)
            rows.append({'variant': variant[self._var_col_name],
                         'gt_fitness': variant[self._fit_col_name],
                         'zs_score': log_zs_score})
        score_df = pd.DataFrame(rows)
        print(score_df['gt_fitness'].corr(score_df['zs_score'], method='spearman'))

        # save score_df
        score_df_path = os.path.join(self.results_dir, f'{os.path.basename(self._input_csv)[:-4]}_score_df.csv')
        score_df.to_csv(score_df_path)
        print(f'ZS score .csv file saved to {score_df_path}.')

    
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
        FlowsiteData(input_csv=lib, **kwargs) 
        
