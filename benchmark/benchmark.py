from mhclovac import predict
from mhclovac.sequence import read_fasta
from mhclovac.utils import load_model
import numpy as np
import re
import matplotlib.pyplot as plt
from matplotlib import rcParams
from sklearn.metrics import auc

rcParams['font.family'] = 'sans-serif'
# rcParams['font.sans-serif'] = ['Tahoma']
rcParams['font.size'] = 12
props = dict(boxstyle='round', facecolor='wheat', alpha=0.6)


bench_fasta = 'data/CD8_epitopes.fsa'

scores = {
    'binding': {
        'positives': [],
        'negatives': []
    },
    'epitope': {
        'positives': [],
        'negatives': []
    },
    'combined': {
        'positives': [],
        'negatives': []
    }
}

random_selection = np.random.randint(low=0, high=1660, size=200)

for i, (sequence_name, sequence) in enumerate(read_fasta(bench_fasta)):

    if i not in random_selection:
        continue

    # True positive, extracted from sequence name
    true_epitope = sequence_name.split()[0]

    # Convert MHC allele to compatible format
    # Example: HLA-A02:01 to HLA-A*02:01 (asterisk is missing)
    allele = sequence_name.split()[1]
    mhc_class = re.findall('HLA-[A-Z]', allele)[0]
    mhc_sub = allele.split(mhc_class)[1]
    mhc_allele = '*'.join([mhc_class, mhc_sub])

    try:
        bmodel, emodel = load_model(mhc_allele)
        if not all([bmodel, emodel]):
            print(mhc_allele, 'skipped')
            continue
    except Exception as e:
        print(e)
        continue

    try:
        predictions = predict(sequence, len(true_epitope), mhc_allele, 'benchmark')
    except Exception as e:
        print(e)
        continue

    for i, row in predictions.iterrows():
        if row['sequence'] == true_epitope:
            scores['binding']['positives'].append(row['binding_score'])
            scores['epitope']['positives'].append(row['epitope_score'])
            scores['combined']['positives'].append(row['combined_score'])
        else:
            scores['binding']['negatives'].append(row['binding_score'])
            scores['epitope']['negatives'].append(row['epitope_score'])
            scores['combined']['negatives'].append(row['combined_score'])


# TPR - true positive rate, fraction of true positives out of all real positives
# FPR - false positive rate, fraction of false positives out of all real negatives

roc_data = {
    'binding': {
        'tpr_array': [],
        'fpr_array': []
    },
    'epitope': {
        'tpr_array': [],
        'fpr_array': []
    },
    'combined': {
        'tpr_array': [],
        'fpr_array': []
    }
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
auc_binding = round(auc(roc_data['binding']['fpr_array'], roc_data['binding']['tpr_array']), 3)
auc_epitope = round(auc(roc_data['epitope']['fpr_array'], roc_data['epitope']['tpr_array']), 3)
auc_combined = round(auc(roc_data['combined']['fpr_array'], roc_data['combined']['tpr_array']), 3)

ax.set_title('ROC curve')
ax.plot(roc_data['binding']['fpr_array'], roc_data['binding']['tpr_array'])
ax.plot(roc_data['epitope']['fpr_array'], roc_data['epitope']['tpr_array'])
ax.plot(roc_data['combined']['fpr_array'], roc_data['combined']['tpr_array'])

ax.grid()
ax.set_xlabel('False positive rate')
ax.set_ylabel('True positive rate')
ax.legend(labels=(
    f'binding AUC = {auc_binding}',
    f'epitope AUC = {auc_epitope}',
    f'combined AUC = {auc_combined}'
))
plt.savefig('results/ROC.png', bbox_inches='tight')
