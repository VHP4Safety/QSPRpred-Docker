

# spock - **S**tandardize, **P**repare, D**ock**

<img src="https://th.bing.com/th/id/OIG.tKEG33W0qiLIBwU0_Fxn?pid=ImgGn"  width="200">

This package contains a minimalistic python API for HTS. Includes pipelines for:
1. Storage, retrieval
2. Standardization
3. Ligand preparation
4. Docking by using [VinaGPU](https://github.com/DeltaGroupNJUPT/Vina-GPU)*   via a Docker image (not really integrated in the package yet)

## *Work in progress*

### Tutorial 

The **tutorial/** folder contains two notebooks:
- standardize.ipynb - creating a MoleculeStore and standardizing ligands from an external library
- ligprep.ipynb - running Schrodinger ligprep on standardized molecules


### Installation

Package can be installed as follows:

```bash
git clone ...
cd spock
pip install -e .


Testing git capabilities.
```
