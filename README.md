# mhc-profiles
Analysis of modeled physicochemical profiles of MHC ligands. 

* [About](#about)
* [Installation](#installation)
* [Usage](#usage)

### About
MhcLovac uses physicochemical index values of amino acid residues to create theoretical models of mhc ligands. 
Based on these models and experimental data we are able to predict binding affinity of MHC-peptide complex with relative confidence.  
We use 4 different physicochemical properties for modelling:
* hydrophobicity
* hydrogen donor bonds
* positive charge
* negative charge

Information about distribution of each of these properties across the peptide's sequence is what is used by MhcLovac to learn to distinguish good from bad binders. 

### Installation
```
pip install mhclovac
```

### Usage
