# MHCLovac 3

MHC class I binding and epitope prediction based on modeled physicochemical properties of peptides.

### Disclaimer
MHCLovac is a result of my personal interests and some free time. 
Although I spent good amount of time researching this topic I would not qualify this project as a thorough scientific research. 
That said, MHCLovac is not too bad in terms of predictions it makes and I plan to improve it further if possible.

### What's new in version 3?
* Epitope prediction in form of epitope probability score. 
* Binding score is back. The ic50 predictions from version 2 are replaced with binding score. Higher score means stronger binding.
* Prediction is carried out by a collection of regression and classification algorithms.
* MHCLovac can now be used as python package. See [example usage](#example-usage) for more info.

### About

MHCLovac is MHC class I binding and epitope prediction tool. 
It uses physicochemical properties of peptides to predict binding affinity and epitope probability (in form of scores). 
One of the main challenges with MHC binding prediction, which MHCLovac aims to solve, is that target peptides don't have to be uniform in length. 
Some alleles allow peptide lengths to span a wide range: H2-Kb epitopes are known to span 7 - 13 residues in length. 
This poses a challenge when creating numerical feature representation of peptides for prediction algorithms. 
MHCLovac solves this by modeling each peptide into a linear, wave-like representation of its physicochemical properties. 
Modeled array can then be scaled up or down to a fixed length allowing MHCLovac to work with a fixed number of features. 
The downside to this approach is that MHCLovac has to assume that peptides bind in linear conformation which may not always be the case.

MHCLovac is trained on data obtained from two sources: 
dataset used for retraining the IEDB class I binding prediction tools [http://tools.iedb.org/main/datasets/](http://tools.iedb.org/main/datasets/) 
and IEDB database [www.iedb.org](www.iedb.org). 
Training results and a list of supported MHC alleles is available in [training/results](training/results) folder.
Trained models are benchmarked using ROC-AUC method. 
Benchmarking method is covered in [benchmark](benchmark) folder.

![roc/auc](https://gitlab.com/stojanovicbg/mhclovac/-/raw/master/benchmark/results/ROC.png)

### Installation

```
pip install mhclovac
```

### Example usage

As command line tool:
```
mhclovac -f example.fasta -m HLA-B*44:02 -l 11
```

As python package:
```python
from mhclovac import predict
from mhclovac.utils import list_mhc_alleles

alleles = list_mhc_alleles()
# returns list of supported MHC alleles

predictions = predict(sequence=['MEIFIEVFSHF', 'ELTLNMCL'], mhc='HLA-B*44:02')
# returns pandas DataFrame with prediction results

```

### Example output
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

### Donate to support my work
If you like this project and wish to support my work you can do so by donating Bitcoin. 
Any amount donated will be appreciated! 

BTC: `bc1qrg7wku5g35kn0qyay4uwzugfmfqwnvz95g54pj`
