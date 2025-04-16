## About
A repository for our paper titled "Substrate-Aware Zero-Shot Predictors for Non-Native Enzyme Activities"

## Environments
* The main environment is `substrate_aware.yml`
* The `coves.yml`, `esmif.yml`, and `plip.yml` files are conda environment for the COVES, ESM-IF, and PLIP zero-shot calcualtion, respectively.
* Frozen versions of the dependencies can be found under `envs/fronzen/`.

### Vina
* Install vina following the instructions in the documentation [here](https://autodock-vina.readthedocs.io/en/latest/installation.html)
* Create a seperate conda environment for vina
```
conda create --name vina_env python=3.9.7 -y
conda activate vina_env
conda install -c conda-forge numpy openbabel pdbfixer scipy rdkit openmm biopython -y
pip install meeko
conda install -c conda-forge mdanalysis
pip install numpy==1.24.4
```
* Run the zero-shot prediction using the following command
```
~/miniconda3/envs/vina/bin/python -m tests.test_zs-vina
```

### LigandMPNN
* Get code and model parameters
```
git clone https://github.com/dauparas/LigandMPNN.git
cd LigandMPNN
bash get_model_params.sh "./model_params"
```
* Modify the paths in `substrate_aware.zs.ligandmpnn` and/or `tests.test_zs-ligandmpnn` accordingly.

### FlowSite
* Get code and moedel parameters
```
git clone https://github.com/HannesStark/FlowSite.git
cd FlowSite
mkdir model_params
cd model_params
```
* Download the model parameters in a zip file based on the repo from [here](https://drive.google.com/file/d/1QGQ6U3BDlEZ682yv7dLo2wbuulOPArSY/view?usp=sharing).
* Unzip as needed and check that there are subfolders called `lf5t55w4` and `b1ribx1a`
* Modify the path in `substrate_aware.zs.ligandmpnn` and/or `tests.test_zs-ligandmpnn` accordingly.
* Set up the conda environment
```
conda create -n flowsite python=3.10
conda activate flowsite
pip install torch==2.1.0+cu121 torchvision==0.16.0+cu121 torchaudio==2.1.0+cu121 -f https://download.pytorch.org/whl/torch_stable.html
pip install torch-scatter -f https://data.pyg.org/whl/torch-2.1.0+cu121.html
pip install torch-sparse -f https://data.pyg.org/whl/torch-2.1.0+cu121.html
pip install torch-cluster -f https://data.pyg.org/whl/torch-2.1.0+cu121.html
pip install torch-spline-conv -f https://data.pyg.org/whl/torch-2.1.0+cu121.html
pip install torch-geometric
pip install rdkit
pip install pyyaml wandb biopython spyrmsd einops biopandas plotly prody tqdm lightning imageio datasets
pip install e3nn
```
* Troubleshoot
```
pip uninstall numpy
pip install numpy==1.23.5
```

## Datasets
The `data/` folder is organized as follows:
* `lib/` contains CSV files with sequence (AAs), fitness, and selectivity (when applicable) data.
    * Each file is named using the format `{enzyme}_{substrate}.csv`.
    * Note: `ParLQ-a` is sometime abbreviated as `ParLQ`
* `seq/` contains FASTA files for each enzyme sequence.
* `structure/` contains the structural files in PDB format.
    * `ParLQ.pdb` is taken from Yang and Lal et al. (2025).
    * `PfTrpB.pdb` corresponds to chain A of PDB entry `5DW0`.
    * `Rma.pdb` corresponds to PDB entry `3CP5`.
    * `apo/` contains apo (substrate-free) structures in PDB format.
    * `clean/` contains cleaned structures (solvent removed) in PDB format.
    * `docked/` contains docked structures in CIF format, generated using AF3 or Chai.
        * See `docked/README.md` for additional information on the docked structures.
* `evmodel/` contains EVcouplings model files (`.model`) for each enzyme, which are used to calculate EVmutation scores.
* `evmodel_dets/` contains detailed metadata for the corresponding EVcouplings models.

## Preprocessing
* Run `preprocess_all` from `substrate_aware.preprocess` to preprocess the data.

## Zero-shot prediction
* Each zero-shot predictor corresponds to a module in `substrate_aware.zs`.
* `run_all_combzs` in `substrate_aware.zs` can be used to combine all generated zero-shot scores for each dataset.

## Analysis
* The analysis scripts are located in `substrate_aware.analysis`.