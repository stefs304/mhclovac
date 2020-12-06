# MHCLovac

[![Downloads](https://pepy.tech/badge/mhclovac)](https://pepy.tech/project/mhclovac)
[![Downloads](https://pepy.tech/badge/mhclovac/week)](https://pepy.tech/project/mhclovac)

MHC binding prediction based on modeled physicochemical properties of peptides.

> Stojanovic S. MHCLovac 4.0: MHC binding prediction based on modeled physicochemical properties of peptides. https://pypi.org/project/mhclovac/4.0

`bitcoin: bc1qrg7wku5g35kn0qyay4uwzugfmfqwnvz95g54pj`

### Table of content
* [About](#about)
* [Materials and methods](#materials-and-methods)
* [Results](#results)
* [Installation](#installation)
* [Usage](#example-usage)
* [References](#)

### About

MHCLovac uses physicochemical properties of peptides to predict binding affinity.
One of the main challenges with MHC binding prediction, which MHCLovac aims to solve, is that target peptides don't have to be uniform in length. 
Some alleles allow peptide lengths to span a wide range: H2-Kb epitopes are known to span 7 - 13 residues in length. 
This poses a challenge when creating numerical feature representation of peptides for prediction algorithms. 
MHCLovac solves this by modeling each peptide into a linear, wave-like representation of its physicochemical properties. 
Modeled array can then be scaled up or down to a fixed length allowing MHCLovac to work with a fixed number of features. 
The downside to this approach is that MHCLovac has to assume that peptides bind in linear conformation which may not always be the case.

### Materials and methods

MHCLovac is trained on data obtained from two sources: 
dataset used for retraining the IEDB class I binding prediction tools [http://tools.iedb.org/main/datasets/](http://tools.iedb.org/main/datasets/) 
and IEDB database [www.iedb.org](www.iedb.org). 
Training results and a list of supported MHC alleles is available in [training/results](training/results) folder.

![mhclovac-physicochemical-profile-peptide](research/figures/mhclovac-modeling-figure.png)

### Results

Trained models are benchmarked using ROC-AUC method. 
Benchmarking method is covered in [benchmark](benchmark) folder.

![roc/auc](https://gitlab.com/stojanovicbg/mhclovac/-/raw/master/benchmark/results/ROC.png)

### Installation

```
pip install mhclovac
```

### Usage

As command line tool:
```
mhclovac -f example.fasta -m HLA-B*44:02 -l 11
```

As python library:
```python
from mhclovac import predict
from mhclovac.utils import list_mhc_alleles

alleles = list_mhc_alleles()
# returns list of supported MHC alleles

predictions = predict(sequence=['MEIFIEVFSHF', 'ELTLNMCL'], mhc='HLA-B*44:02')
# returns pandas DataFrame with prediction results

```

Example output:
```
 sequence          mhc  peptide_length           sequence_name  binding_score  epitope_score  combined_score
 MEIFIEVFSHF  HLA-B*44:02              11  MEIFIEVFSHF HLA-B44:02       0.523205       0.965484        1.488688
 EIFIEVFSHFL  HLA-B*44:02              11  MEIFIEVFSHF HLA-B44:02       0.087188       0.512132        0.599320
 IFIEVFSHFLL  HLA-B*44:02              11  MEIFIEVFSHF HLA-B44:02       0.039142       0.159362        0.198503
 FIEVFSHFLLQ  HLA-B*44:02              11  MEIFIEVFSHF HLA-B44:02       0.114877       0.264553        0.379430
 IEVFSHFLLQL  HLA-B*44:02              11  MEIFIEVFSHF HLA-B44:02       0.317922       0.964168        1.282090
```

Columns:
1. `sequence` 
2. `sequence_name` - Fasta sequence name or name provided by `-n` argument
3. `peptide_length`
4. `mhc` - MHC allele
5. `binding_score` - Higher score means better binding
6. `epitope_score` - Higher score means a better epitope
7. `combined_score` - Sum of binding and epitope scores if both are available

### References
1. Vita R, Mahajan S, Overton JA, Dhanda SK, Martini S, Cantrell JR, Wheeler DK, Sette A, Peters B. The Immune Epitope Database (IEDB): 2018 update. Nucleic Acids Res. 2018 Oct 24. doi: 10.1093/nar/gky1006. [Epub ahead of print] PubMed PMID: [30357391](https://www.ncbi.nlm.nih.gov/pubmed/30357391).
2. Benchmarking

