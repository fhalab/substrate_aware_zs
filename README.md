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
```
* Run the zero-shot prediction using the following command
```
~/miniconda3/envs/vina/bin/python -m tests.test_zs-vina
```

### LigandMPNN

(For developers instructions: [https://github.com/dauparas/LigandMPNN](https://github.com/dauparas/LigandMPNN))

```python
git clone https://github.com/dauparas/LigandMPNN.git
cd LigandMPNN
bash get_model_params.sh "./model_params"

conda create -n ligandmpnn_env python=3.11.8
pip3 install -r requirements.txt
```

- install PyTorch according to your build ([https://pytorch.org/get-started/locally/](https://pytorch.org/get-started/locally/))
