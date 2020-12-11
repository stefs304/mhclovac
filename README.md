# MHCLovac :serbia:

<!--
  Title: MHCLovac
  Description: MHC binding prediction based on modeled physicochemical properties of peptides
  Author: Stefan Stojanovic
  Keywords: mhc, binding, predcition, ligand, immuno, physicochemical, peptides, modeling
  -->

[![Downloads](https://pepy.tech/badge/mhclovac)](https://pepy.tech/project/mhclovac)
[![Downloads](https://pepy.tech/badge/mhclovac/week)](https://pepy.tech/project/mhclovac)

MHC binding prediction based on modeled physicochemical properties of peptides.

`bitcoin: bc1qrg7wku5g35kn0qyay4uwzugfmfqwnvz95g54pj`

### What's new in version 4
* Training data is better curated, duplicate samples are removed and MHC allele names are standardized.
* Models are trained on MHC binding affinity data and eluted ligand data.
* Allele specific set of physicochemical properties, determined based on individual prediction ability.
* Ligand prediction from the previous version was removed due to training data not being properly curated. 
* `binding_score` is reported as the log transformed binding affinity: `1 - log50k(ic50)`, as per NetMHC standard.


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
 peptide          mhc  peptide_length           sequence_name  binding_score
 MEIFIEVFSHF  HLA-B*44:02              11  MEIFIEVFSHF HLA-B44:02       0.523205
 EIFIEVFSHFL  HLA-B*44:02              11  MEIFIEVFSHF HLA-B44:02       0.087188
 IFIEVFSHFLL  HLA-B*44:02              11  MEIFIEVFSHF HLA-B44:02       0.039142
 FIEVFSHFLLQ  HLA-B*44:02              11  MEIFIEVFSHF HLA-B44:02       0.114877
 IEVFSHFLLQL  HLA-B*44:02              11  MEIFIEVFSHF HLA-B44:02       0.317922
```

Columns:
1. `peptide` 
2. `sequence_name` - Fasta sequence name or name provided by `-n` argument
3. `peptide_length` -
4. `mhc` - MHC allele name
5. `binding_score` - Higher score means better binding


