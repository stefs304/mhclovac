from mhclovac.sequence import read_fasta, chop_sequence
import mhclovac
import pandas as pd
import mhcnames
import joblib
import numpy as np


N_SAMPLES = 200
RANDOM_SEED = 0
bench_fasta = '../data/CD8_epitopes.fsa'


def get_mhclovac_frank(predictions, epitope):
    count = 0.0
    for i, row in predictions.iterrows():
        if row['peptide'] == epitope:
            break
        count += 1
    return count / len(predictions)


scores = {
    'true_positive_score': [],
    'true_negative_score': []
}
frank_scores = []

np.random.seed(RANDOM_SEED)
random_selection = np.random.randint(0, 1659, N_SAMPLES)

for i, (sequence_name, sequence) in enumerate(read_fasta(bench_fasta)):
    if i not in random_selection:
        continue

    print(i, sequence_name)
    true_epitope = sequence_name.split()[0]
    mhc_allele = mhcnames.normalize_allele_name(sequence_name.split()[1])

    try:
        peptides = chop_sequence(sequence, len(true_epitope))
        predictions = mhclovac.predict(peptides=peptides, mhc_allele=mhc_allele)
        frank = get_mhclovac_frank(predictions, true_epitope)
        frank_scores.append(frank)
        for _, row in predictions.iterrows():
            if row['peptide'] == true_epitope:
                scores['true_positive_score'].append(row['binding_score'])
            else:
                scores['true_negative_score'].append(row['binding_score'])
    except Exception as e:
        print('mhclovac', e)

data = pd.DataFrame()
data['mhclovac'] = frank_scores
data.to_csv('./results/benchmarking_frank.csv', index=False)

roc_data = {
    'tpr_array': [],
    'fpr_array': []
}

t_min = min([min(scores['true_positive_score']), min(scores['true_negative_score'])])
t_max = max([max(scores['true_positive_score']), max(scores['true_negative_score'])])

roc_data['tpr_array'].append(1.0)
roc_data['fpr_array'].append(1.0)

for t in sorted(scores['true_positive_score']):
    true_pos = len([x for x in scores['true_positive_score'] if x >= t])
    false_pos = len([x for x in scores['true_negative_score'] if x >= t])
    tpr = float(true_pos) / len(scores['true_positive_score'])
    fpr = float(false_pos) / len(scores['true_negative_score'])
    roc_data['tpr_array'].append(tpr)
    roc_data['fpr_array'].append(fpr)

roc_data['tpr_array'].append(0.0)
roc_data['fpr_array'].append(0.0)


joblib.dump(roc_data, 'results/benchmarking_roc.gzip', compress=('gzip', 5))

