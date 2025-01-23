import subprocess
import os
import re
import ipdb

from glob import glob
from copy import deepcopy
from typing import Union

from REVIVAL.preprocess import ZSData


class PyrosettaData(ZSData):


    def __init__(
        self,
        input_csv: str,
        scale_fit: str = "parent",
        chain_id: str = "A",
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        mut_col_name: str = "mut",
        pos_col_name: str = "pos",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        structure_dir: str = "data/structure",
        preprocessing_dir: str ="data/tmp",
        preprocessed_dir: str ="data/preprocessed",
        rosetta_param_script_path: str = "./rosetta/source/scripts/python/public/generic_potential/mol2genparams.py"
    ):
        super().__init__(            
            input_csv=input_csv,
            scale_fit=scale_fit,
            combo_col_name=combo_col_name,
            var_col_name=var_col_name,
            mut_col_name=mut_col_name,
            pos_col_name=mut_col_name,
            seq_col_name=seq_col_name,
            fit_col_name=fit_col_name,
            structure_dir=structure_dir, 
            )
        
        self.chain_id = chain_id
        self.rosetta_param_script_path = rosetta_param_script_path
        self.structure_dir = structure_dir
        self.preprocessing_dir = preprocessing_dir
        self.preprocessed_dir = preprocessed_dir

        self.preprocess_inputs()


    def split_enzyme_and_ligand(self, pdb_lines, struc_path):
        ligand_lines = []
        enzyme_lines = []
        for pdb_line in pdb_lines:
            if (pdb_line.startswith('ATOM') or  pdb_line.startswith('HETATM'))  and pdb_line[17:20].strip()=='LIG':
                ligand_lines.append(pdb_line)
            elif (pdb_line.startswith('ATOM') or  pdb_line.startswith('HETATM'))  and not pdb_line[17:20].strip()=='LIG':
                enzyme_lines.append(pdb_line)

        ligand_pdb_file = os.path.basename(struc_path)[:-4] + '_ligand' + '.pdb'
        ligand_pdb_filepath = os.path.join(self.preprocessing_dir, ligand_pdb_file)

        enzyme_pdb_file = os.path.basename(struc_path)[:-4] + '_enzyme' + '.pdb'
        enzyme_pdb_filepath = os.path.join(self.preprocessing_dir, enzyme_pdb_file)

        with open(ligand_pdb_filepath, 'w') as ligpdb:
            ligpdb.writelines(ligand_lines)
            print(f'INFO: Ligand atoms succesfully extracted and saved as .pdb file at {ligand_pdb_file}.')
        
        with open(enzyme_pdb_filepath, 'w') as enzpdb:
            enzpdb.writelines(enzyme_lines)
            print(f'INFO: Enzyme atoms succesfully extracted and saved as .pdb file at {enzyme_pdb_file}.')
        
        return enzyme_pdb_filepath, ligand_pdb_filepath
    
    
    def pdb2mol2(self, ligand_pdb_path):
        out_path = os.path.join(self.preprocessing_dir, os.path.basename(ligand_pdb_path))[:-4] + '.mol2'
        cli = ['obabel',
               ligand_pdb_path,
               '-O',
               out_path]
        cli_return = subprocess.run(cli, capture_output=True, text=True)
        print(f'INFO: Succesfully generated .mol2 file at {out_path}')
        return out_path


    def generate_charges(self, ligand_mol2_path):

        # if multiple units are present, they must be split in separate mol2 for antechamber
        cli = ['obabel',
                ligand_mol2_path,
                '-O',
                ligand_mol2_path[:-5] + '_unit.mol2',
                '-m',
                '--separate'       
            ]
        cli_return = subprocess.run(cli, capture_output=True, text=True)
        print(f'INFO: Succesfully ligand mol2 file along its units.')

        # generate the am1bcc charges for both units separately using antechamber
        tmp_files_list = os.listdir(self.preprocessing_dir)
        ligand_name = os.path.basename(ligand_mol2_path)
        unit_list = [unit for unit in tmp_files_list if unit.startswith(ligand_name[:-5]+'_unit') and not unit.endswith('_am1bcc.mol2') and not unit.endswith('.params')]
        charge_list = [0, -3] #TODO: set up system for passing this
        charged_unit_paths = []
        for i, unit_name in enumerate(unit_list):
            unit_path = os.path.join(self.preprocessing_dir, unit_name)
            out_path = os.path.join(self.preprocessing_dir, os.path.basename(unit_path))[:-5] + '_am1bcc.mol2'        
            cli = ['antechamber',
                '-i',
                unit_path,
                '-fi',
                'mol2',
                '-o',
                out_path,
                '-fo', 
                'mol2',
                '-c', 
                'bcc',       # defines chargemethod bcc = am1bcc
                '-nc',       # defines net charge
                str(charge_list[i])
            ]

            cli_return = subprocess.run(cli, capture_output=True, text=True)
            charged_unit_paths.append(out_path)
            print(f'INFO: Sucessfully generated am1bcc charges for {out_path}.')

        # merging the two units back together #
        cli = ['obabel'] + charged_unit_paths + [
               '-O',
               unit_path[:-11]+'_am1bcc.mol2', '--combine', 'm']
        cli_return = subprocess.run(cli, capture_output=True, text=True)
        print(f'INFO: Succesfully recombined units into a am1bcc_charged.mol2 file.')
        
        return unit_list


    def generate_params(self, unit_list):
        for unit_name in unit_list: 
            unit_am1bcc_path = os.path.join(self.preprocessing_dir, unit_name)
            cli = ['python',
                self.rosetta_param_script_path,
                '-s',
                unit_am1bcc_path,
                '--outdir',
                self.preprocessing_dir,
                '--no_pdb'
            ]
            cli_return = subprocess.run(cli, capture_output=True, text=True)


    def ligands_into_enzyme(self, unit_list, enzyme_pdb_path):
        """
        Only works in unix-based systems (i.e. Linux/Ubuntu or macOS)
        """
        merge_list = [os.path.join(self.preprocessing_dir, unit) for unit in unit_list]
        cli = f"cat {enzyme_pdb_path} {' '.join(merge_list)} > {enzyme_pdb_path[:-4]}_merged.pdb"
        cli_return = subprocess.run(cli, capture_output=True, text=True, shell=True)
        ipdb.set_trace()


    def synchronize_pdb_with_params(self, pdb_lines, params_path_list, pdb_output_path):
        chain_names = ['X', 'Y', 'Z']  # Alter if more than three ligands are needed

        # find the last protein_line_index and last_residue_number
        for i, pdb_line in enumerate(pdb_lines):
            if pdb_line.startswith('TER'):
                last_protein_line = pdb_lines[i-1]
                last_protein_line_index = int(last_protein_line[7:11])
                last_residue_number = int(last_protein_line[23:26])
                break

        # pre-process pdb atom numbering if necessary
        elements_list = []
        for pdb_line in pdb_lines:
            if (pdb_line.startswith('ATOM') or  pdb_line.startswith('HETATM'))  and pdb_line[17:20].strip()=='LIG':
                elements_list.append(pdb_line[11:16].strip())
        elements_set = set(elements_list)

        # initialize counter dict
        element_counter_dict = {}
        for element in elements_set:
            element_counter_dict[element] = 0

        element_max_dict = {}
        for element in elements_set:
            element_max_dict[element] = 0
        for i, pdb_line in enumerate(pdb_lines):
            if (pdb_line.startswith('ATOM') or  pdb_line.startswith('HETATM'))  and pdb_line[17:20].strip()=='LIG':
                element = pdb_line[11:16].strip()
                element_max_dict[element] += 1

        for i, pdb_line in enumerate(pdb_lines):
            if (pdb_line.startswith('ATOM') or  pdb_line.startswith('HETATM'))  and pdb_line[17:20].strip()=='LIG':
                element = pdb_line[11:16].strip()
                if element_max_dict[element] != 1:
                    element_counter_dict[element] += 1
                    new_element = element + str(element_counter_dict[element])
                else:
                    new_element = element
                pdb_linelist = pdb_line.split()
                pdb_line = "{:<6}{:>5}  {:<4}{:>3} {:>1}{:>4}{:>11}{:>8}{:>8} {:>6}{:>6}{:>12}\n".format(
                                                        pdb_linelist[0],  # Record type
                                                        pdb_linelist[1],  # Atom serial number
                                                        new_element,      # Atom name
                                                        pdb_linelist[3],  # Residue name
                                                        pdb_linelist[4],  # Chain identifier
                                                        pdb_linelist[5],  # Residue sequence number
                                                        pdb_linelist[6],  # X coordinate
                                                        pdb_linelist[7],  # Y coordinate
                                                        pdb_linelist[8],  # Z coordinate
                                                        pdb_linelist[9],  # Occupancy
                                                        pdb_linelist[10], # Temperature factor
                                                        pdb_linelist[11]  # Element symbol
                                                        )
                pdb_lines[i] = pdb_line
            
        # Iterate over lines in the parameter files
        pdb_ligand_lines = []
        for i, path in enumerate(params_path_list):
            with open(path, 'r') as param:
                ligand_name = None
                
                for lig_line in param:
                    if lig_line.startswith('NAME'):
                        ligand_name = lig_line.split()[1]
                        break
                if ligand_name == None:
                    print(f'ERROR: Could not find ligand name for {path}.')
                
                for lig_line in param:
                    if lig_line.startswith('ATOM'):
                        lig_line_list = lig_line.split()
                        lig_atom_name = lig_line_list[1]

                                # obtain the last line before the first ter


                        # find the line in the pdb that corresponds to this ligand atom and correct it
                        found_flag = False
                        for pdb_line in pdb_lines:
                            
                            if (pdb_line.startswith('ATOM') or  pdb_line.startswith('HETATM')) and (pdb_line[13:16].strip()==lig_atom_name or pdb_line[13:16].strip()==lig_atom_name[0]+'0'+lig_atom_name[-1]) and pdb_line[17:20].strip()=='LIG':
                                
                                found_flag = True
                                pdb_linelist = pdb_line.split()
                                if len(pdb_linelist) != 12:
                                    ipdb.set_trace()    
                                new_pdb_line = "{:<6}{:>5}  {:<4}{:>3} {:>1}{:>4}{:>11}{:>8}{:>8} {:>6}{:>6}{:>12}\n".format(
                                                    pdb_linelist[0],  # Record type
                                                    last_protein_line_index + len(pdb_ligand_lines) + 1,  # Atom serial number
                                                    lig_atom_name,  # Atom name
                                                    ligand_name,  # Residue name
                                                    pdb_linelist[4],  # Chain identifier
                                                    last_residue_number + i + 1,  # Residue sequence number
                                                    pdb_linelist[6],  # X coordinate
                                                    pdb_linelist[7],  # Y coordinate
                                                    pdb_linelist[8],  # Z coordinate
                                                    pdb_linelist[9],  # Occupancy
                                                    pdb_linelist[10], # Temperature factor
                                                    pdb_linelist[11]  # Element symbol
                                                    )
                                pdb_ligand_lines.append(new_pdb_line)

                        if found_flag == False:
                            print(f'Failed to find the atom {lig_atom_name} of ligand {ligand_name} in the PDB file.')

                
                # find the end of the residue atoms
                pdb_lines_to_keep = []
                for pdb_line in pdb_lines:
                    pdb_lines_to_keep.append(pdb_line)
                    if pdb_line.startswith('TER'):
                        break

                # add the lines that contain the ligand atoms
                for pdb_ligand_line in pdb_ligand_lines:
                    pdb_lines_to_keep.append(pdb_ligand_line)
                pdb_lines_to_keep.append('TER\n')
                
                with open(pdb_output_path, 'w') as output_file:
                    output_file.writelines(pdb_lines_to_keep)
                print(f'SUCESS: Saved reformated PDB under {pdb_output_path}')

 


        
    def preprocess_inputs(self):
        """
        Input .pdb file is expected to hold an enzyme with docked ligand.
        Ligand requires correct bond orders and hydrogens to be defined. 
        The active site must not contain single atoms (e.g. ions), or molecules with less than four atoms.
        """
        for struc_filename in  os.listdir(self.structure_dir):
            struc_path = os.path.join(self.structure_dir, struc_filename)
            with open(struc_path, 'r') as PDB:
                pdb_lines = PDB.readlines()
            
            # extract the ligand and save as ligand.pdb
            enzyme_pdb_path, ligand_pdb_path = self.split_enzyme_and_ligand(pdb_lines, struc_path)

            # convert ligand to mol2
            ligand_mol2_path = self.pdb2mol2(ligand_pdb_path)

            # generate am1bcc charges for ligand.mol from antechamber
            unit_list = self.generate_charges(ligand_mol2_path)

            # generate params file from ligand.mol2
            params_paths = self.generate_params(unit_list)

            # (opt.: convert the ligand.mol2 back to a pdb)

            # merge ligand.pdb with protein.pdb
            merged_pdb = self.ligands_into_enzyme(unit_list, enzyme_pdb_path)

            # synchronize ligand.params and ligand.pdb
            synchronized_pdb = self.synchronize_pdb_with_params(pdb_lines, params_paths, self.preprocessed_dir)

            # Run PyRosetta GALigandDock
            self.inference()
    
    def inference():
        pass



def run_pyrosetta_pipeline(pattern: Union[str, list] = None, kwargs: dict = {}):
    if isinstance(pattern, str):
        lib_list = glob(pattern)
    else:
        lib_list = deepcopy(pattern)

    for lib in lib_list:
        print(f"Running Pyrosetta Pipeline for {lib}...")
        PyrosettaData(input_csv=lib, **kwargs)

dummy_pattern = ['_.csv']
dummy_kwargs = {"structure_dir": "/disk2/lukas/dummy_pyrosetts_data/input",
                "preprocessing_dir": "/disk2/lukas/dummy_pyrosetts_data/tmp",
                "preprocessed_dir": "/disk2/lukas/dummy_pyrosetts_data/preprocessed",
                "rosetta_param_script_path": "/disk2/lukas/rosetta/source/scripts/python/public/generic_potential/mol2genparams.py"}

run_pyrosetta_pipeline(pattern=dummy_pattern, kwargs=dummy_kwargs)