# REVIVAL2
A repository for the REVIVAL2 project, applying SSMuLA to non-native functions.

## Environments
* The main environment is `REVIVAL.yml`
* The `coves.yml` and `esmif.yml` files are conda environment for the COVES and ESMIF zero-shot calcualtion, respectively.

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
* Modify the paths in `REVIVAL.zs.ligandmpnn` and/or `tests.test_zs-ligandmpnn` accordingly.

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
* Modify the path in `REVIVAL.zs.ligandmpnn` and/or `tests.test_zs-ligandmpnn` accordingly.
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
pip install pyyaml wandb biopython spyrmsd einops biopandas plotly prody tqdm lightning imageio
pip install e3nn
```
* Troubleshoot
```
pip uninstall numpy
pip install numpy==1.23.5
```