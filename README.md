# MHCLovac

<!--
  Title: MHCLovac
  Description: MHC binding prediction based on modeled physicochemical properties of peptides
  Author: Stefan Stojanovic
  Keywords: mhc, binding, predcition, ligand, immuno, physicochemical, peptides, modeling
  -->

[![Downloads](https://pepy.tech/badge/mhclovac)](https://pepy.tech/project/mhclovac)
[![Downloads](https://pepy.tech/badge/mhclovac/week)](https://pepy.tech/project/mhclovac)

MHC binding prediction based on modeled physicochemical properties of peptides.

> Stojanovic S. MHCLovac 4.0: MHC binding prediction based on modeled physicochemical properties of peptides. https://pypi.org/project/mhclovac/4.0

`bitcoin: bc1qrg7wku5g35kn0qyay4uwzugfmfqwnvz95g54pj`

### Table of content
* [Introduction](#introduction)
* [Materials and methods](#materials-and-methods)
  * [Modeling physicochemical properties](#modeling-physicochemical-properties)
* [Results](#results)
* [Installation](#installation)
* [Usage](#example-usage)
* [References](#references)

### Introduction

MHCLovac is a binding prediction tool which implements the newly developed method for peptide representation.
This method is aimed at overcoming problems caused by variable length of MHC ligands.
To achieve this goal I make a distinction between physicochemical properties of underlying amino acid residues and those of a peptide as a whole, under the following assumption: 
physicochemical properties of adjacent residues in a peptide have additive effect on the local properties of the peptide as a whole and properties of a single residue affect the properties of the peptide at neighboring positions. 
This assumption provides basis for obtaining physicochemical profiles of peptides, which allows for their direct comparison regardless of differences in length. 
The prediction accuracy of MHCLovac comparable to that of the state-of-the-art prediction tools and its ability to make predictions for peptides of arbitrary lengths exemplify the practical advantage of this method. 

### Materials and methods

To train MHCLovac I used two datasets obtained from IEDB<sup>[**[1]**](#references)</sup>. 
The first is a dataset used for retraining the IEDB class I binding prediction tools, which contains only the binding affinity measurements ([http://tools.iedb.org/main/datasets](http://tools.iedb.org/main/datasets/)). 
The second is exported MHC ligand dataset comprised of eluted MHC ligands ([https://www.iedb.org/database_export_v3.php](https://www.iedb.org/database_export_v3.php)). 
This dataset contains both quantitative binding affinity measurements and qualitative measurements. 
I combined data from two dataset where quantitative measurements were available into one hybrid dataset used for training the binding prediction model.
For ligand prediction models I used the second dataset. 
The list of all available physicochemical properties and the corresponding amino acid index data was obtained from Aaindex<sup>[**[2]**](#references)</sup> database ([https://www.genome.jp/aaindex](https://www.genome.jp/aaindex)). 

#### Modeling physicochemical properties
The peptide of length `L` is model by creating a vector `S` containing `L*m + 2*m` data points, where `m` is an arbitrary multiplier. 
Each amino acid residue gets a designated slice of the vector `S` corresponding to its relative position in the sequence. 
The i<sup>th</sup> amino acid residue A<sub>i</sub> is modeled as a Gaussian curve G(A<sub>i</sub>) scaled by the corresponding index value G(A<sub>i</sub>) * I<sub>A</sub> (Figure 1.a, dashed colored lines). 
Each G(A<sub>i</sub>) is spans one neighboring slice on each side and is shaped by sigma parameter with default value of 0.8.
The physicochemical profile of peptide is obtained by taking the sum of individually modeled residues (Figure 1.a, black solid line). 
The leading and trailing slices, corresponding to `+ 2*m` term in the first expression, serve the role of placeholders when modeling the first and the last residue. 
These slices are optionally removed to produce the final vector of length `L*m` (not shown in the figure).

![mhclovac-physicochemical-profile-peptide](research/figures/mhclovac-modeling-figure.png)

Prior to modeling each physicochemical property index is normalized to range [-1, 1].
This is done in order for values to reflect relations between different residues on opposite sides of the spectrum. 

Once the profile is obtained it can be further reduced to fixed number of discrete values. 
To give an example, I model two ligands of HLA-A*02:01, one an 8-mer and the other an 11-mer (Figure 1.b and 1.c).
Their profiles are reduced to 10 discrete points each, by sampling the modeled profiles at equal intervals (Figure 1.d and 1.e). 
This way the peptides of different lengths can be reduced to the same fixed number of features. 
>From the Figure 1 it is evident that despite different sequences and lengths, peptides LLDVTAAV and FLFDGSPTYVL have very similar discrete hydrophobicity profiles which can be directly compared. 



#### Prediction models

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
2. Kawashima, S., Pokarowski, P., Pokarowska, M., Kolinski, A., Katayama, T., and Kanehisa, M.; AAindex: amino acid index database, progress report 2008. Nucleic Acids Res. 36, D202-D205 (2008). PMID:[17998252](https://www.ncbi.nlm.nih.gov/sites/entrez?Db=pubmed&Cmd=ShowDetailView&TermToSearch=17998252&ordinalpos=9&itool=EntrezSystem2.PEntrez.Pubmed.Pubmed_ResultsPanel.Pubmed_RVDocSum)

