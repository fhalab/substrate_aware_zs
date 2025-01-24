# Environments Set-Up


### FlowSite

(For developers instructions: [https://github.com/HannesStark/FlowSite](https://github.com/HannesStark/FlowSite))

```python
conda create -c conda-forge -n flowsite rdkit python==3.11.8
pip install torch torchvision torchaudio
pip install torch_geometric
pip install pyyaml wandb biopython spyrmsd einops biopandas plotly prody tqdm lightning imageio
pip install e3nn
pip install torch_scatter torch_sparse torch_cluster -f [https://data.pyg.org/whl/torch-2.1.0+cu121.html](https://data.pyg.org/whl/torch-2.1.0+cu121.html)
```

### Hydro

First install the plip package for active site extraction (for developers instructions: [https://github.com/pharmai/plip/tree/master](https://github.com/pharmai/plip/tree/master))

```python
conda create -n hydro_env python==3.11.11
conda install -c conda-forge openbabel=3.1.1

git clone https://github.com/pharmai/plip.git
cd ./plip
python setup.py install
```

If plip does throw errors upon execution, try commenting out the body of the function int32_to_negative() in plip/plip/baseic/supplemental.py (the function should simply return its input)

The remaining packages to run the [hydro.py](http://hydro.py) script:

```python
conda install -c conda-forge numpy==2.2.0
conda install -c conda-forge pandas==2.2.3
conda install -c conda-forge asttokens==3.0.0
conda install -c conda-forge xmltodict==0.14.2
conda install -c conda-forge tqdm==4.67.1
conda install -c conda-forge scipy==1.14.1
conda install -c conda-forge ipdb==0.13.13
conda install -c conda-forge rdkit==2024.09.3
conda install -c conda-forge openbabel==3.1.1
```