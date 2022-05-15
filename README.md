# `molviewer`

*A Python package for molecular data retrieval and visualisation.*

[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=3600&style=flat-square)](https://github.com/ashtonjesse/molviewer)
[![Build](https://github.com/ashtonjesse/molviewer/actions/workflows/python-package.yml/badge.svg)](https://github.com/ashtonjesse/molviewer/actions)
[![codecov](https://codecov.io/gh/ashtonjesse/molviewer/branch/master/graph/badge.svg?token=DP9RT8XBMI)](https://codecov.io/gh/ashtonjesse/molviewer)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 🔩 Core module

The base logic is handled by the [`molecules`](https://github.com/ashtonjesse/molviewer/blob/master/src/molviewer/molecules.py)
module. It implements two classes `Macromolecule` and `ChemicalMolecule`, 
which retrieve molecular data from the [rcsb](https://www.rcsb.org) and 
[ChEMBL](https://www.ebi.ac.uk/chembl/) databases, respectively. The 
package is intended to be used in a Jupyter Notebook. Users can create class
instances by specifying either a protein data bank id (for `Macromolecule`) 
or a chembl id (for `ChemicalMolecule`).
Both 
classes expose a `show` method for displaying 3D molecule models in a 
Jupyter Notebook and a `save` method for saving molecule data to file. 
`Macromolecule` models are generated using methods available in the 
[`nglview`](https://github.com/nglviewer/nglview) Python package, 
whereas 3D conformers are generated for the `ChemicalMolecule` class using 
methods available in the [`rdkit`](https://www.rdkit.org/) Python package. 

## ⚡ Installation instructions

Create a virtual environment and activate it using the following command in 
the terminal:
```python
.\Scripts\activate
```
If loading scripts is disabled then first run this command in the terminal:
```python
Set-ExecutionPolicy Unrestricted -Scope Process
```
Then install the `molviewer` package:
```python
py -m pip install -i https://test.pypi.org/pypi/ --extra-index-url 
https://pypi.org/simple molviewer==0.0.12
```
Run `Jupyter` with the following command:
```python
jupyter notebook
```
Create a new `notebook` and import the `molviewer` package:
```python
from molviewer.molecules import ChemicalMolecule, Macromolecule
```
    
## 💾 A note on saving .pdb files
While the PDBList class of the Biopython package can retrieve .pdb files from the www.rcsb.org ftp server, these files are retrieved with the '.ent' extension. Instead, we opted to retrieve .pdb files directly from 'https://files.rcsb.org/download/' as these have the correct '.pdb' extension. Hence the final release of this package does not include Biopython in its dependencies.  
