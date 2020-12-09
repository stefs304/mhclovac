from mhclovac import predict
from mhclovac.sequence import read_fasta, chop_sequence
from mhclovac.utils import load_model
import numpy as np
import re
import matplotlib.pyplot as plt
from matplotlib import rcParams
from sklearn.metrics import auc
from mhcflurry import Class1PresentationPredictor
predictor = Class1PresentationPredictor.load()


rcParams['font.family'] = 'sans-serif'
# rcParams['font.sans-serif'] = ['Tahoma']
rcParams['font.size'] = 12
props = dict(boxstyle='round', facecolor='wheat', alpha=0.6)


bench_fasta = '../data/CD8_epitopes.fsa'

scores = {
    'mhclovac': {
        'positives': [],
        'negatives': []
    },
    'mhcflurry': {
        'positives': [],
        'negatives': []
    }
}

random_selection = np.random.randint(low=0, high=1660, size=200)

for i, (sequence_name, sequence) in enumerate(read_fasta(bench_fasta)):

    if i not in random_selection:
        continue

    print(sequence_name)
    # True positive, extracted from sequence name
    true_epitope = sequence_name.split()[0]

    # Convert MHC allele to compatible format
    # Example: HLA-A02:01 to HLA-A*02:01 (asterisk is missing)
    allele = sequence_name.split()[1]
    mhc_class = re.findall('HLA-[A-Z]', allele)[0]
    mhc_sub = allele.split(mhc_class)[1]
    mhc_allele = '*'.join([mhc_class, mhc_sub])

    try:
        model = load_model(mhc_allele)['binding_model']
        if not model:
            print(mhc_allele, 'skipped')
            continue
    except Exception as e:
        print(e)
        continue

    try:
        peptides = chop_sequence(sequence=sequence, peptide_length=len(true_epitope))
        mhclovac_predictions = predict(peptides,  mhc_allele)
        mhcflurry_predictions = predictor.predict(peptides=peptides, alleles=[allele])
    except Exception as e:
        print(e)
        continue

    for i, row in mhclovac_predictions.iterrows():
        if row['sequence'] == true_epitope:
            scores['mhclovac']['positives'].append(row['binding_score'])
        else:
            scores['mhclovac']['negatives'].append(row['binding_score'])

    for i, row in mhcflurry_predictions.iterrows():
        if row['peptide'] == true_epitope:
            scores['mhcflurry']['positives'].append(row['presentation_score'])
        else:
            scores['mhcflurry']['negatives'].append(row['presentation_score'])

# TPR - true positive rate, fraction of true positives out of all real positives
# FPR - false positive rate, fraction of false positives out of all real negatives

roc_data = {
    'mhclovac': {
        'tpr_array': [],
        'fpr_array': []
    },
    'mhcflurry': {
        'tpr_array': [],
        'fpr_array': []
    },
}

for key in scores:
    t_min = min([min(scores[key]['positives']), min(scores[key]['negatives'])])
    t_max = max([max(scores[key]['positives']), max(scores[key]['negatives'])])

    roc_data[key]['tpr_array'].append(1.0)
    roc_data[key]['fpr_array'].append(1.0)

    for t in sorted(scores[key]['positives']):
        true_pos = len([x for x in scores[key]['positives'] if x >= t])
        false_pos = len([x for x in scores[key]['negatives'] if x >= t])
        tpr = float(true_pos) / len(scores[key]['positives'])
        fpr = float(false_pos) / len(scores[key]['negatives'])
        roc_data[key]['tpr_array'].append(tpr)
        roc_data[key]['fpr_array'].append(fpr)

    roc_data[key]['tpr_array'].append(0.0)
    roc_data[key]['fpr_array'].append(0.0)

fig, ax = plt.subplots(1,1, figsize=(5, 5))
auc_mhclovac = round(auc(roc_data['mhclovac']['fpr_array'], roc_data['mhclovac']['tpr_array']), 3)
auc_mhcfurry = round(auc(roc_data['mhcflurry']['fpr_array'], roc_data['mhcflurry']['tpr_array']), 3)

ax.set_title('ROC curve')
ax.plot(roc_data['mhclovac']['fpr_array'], roc_data['mhclovac']['tpr_array'])
ax.plot(roc_data['mhcflurry']['fpr_array'], roc_data['mhcflurry']['tpr_array'])

ax.grid()
ax.set_xlabel('False positive rate')
ax.set_ylabel('True positive rate')
ax.legend(labels=(
    f'MHCLovac AUC = {auc_mhclovac}',
    f'MHCflurry AUC = {auc_mhcfurry}'
))

plt.savefig('figures/binding-roc-auc.png', bbox_inches='tight')
