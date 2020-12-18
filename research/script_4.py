import subprocess
from mhclovac.sequence import read_fasta, chop_sequence
import mhclovac
import pandas as pd
import mhcnames
import os
from io import StringIO
import numpy as np


N_SAMPLES = 200
RANDOM_SEED = 0
IEDB_SCRIPT_PATH = '/home/stefan/mhc_i/src/predict_binding.py'
bench_fasta = '../data/CD8_epitopes.fsa'
IEDB_METHODS = [
    'ann',
    'comblib_sidney2008',
    'consensus',
    # 'IEDB_recommended',
    'netmhcpan_ba',
    'netmhcpan_el',
    'smm',
    'smmpmbec',
    'pickpocket',
    'netmhccons',
    'netmhcstabpan'
]


def get_iedb_method_frank(fasta, mhc, epitope, method):
    peptide_length = str(len(epitope))
    first = subprocess.Popen(['/bin/echo', '-e', fasta], stdout=subprocess.PIPE)
    second = subprocess.Popen(['python', IEDB_SCRIPT_PATH, method, mhc, peptide_length], stdin=first.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = second.communicate()
    dataframe = pd.read_csv(StringIO(str(stdout, 'utf-8')), sep='\t')
    count = 0.0
    for i, row in dataframe.iterrows():
        if row['peptide'] == epitope:
            break
        count += 1
    return count / len(dataframe)


def get_mhclovac_frank(sequence, mhc, epitope):
    peptides = chop_sequence(sequence, len(epitope))
    dataframe = mhclovac.predict(sequence=peptides, mhc=mhc)
    count = 0.0
    for i, row in dataframe.iterrows():
        if row['peptide'] == epitope:
            break
        count += 1
    return count / len(dataframe)


frank_scores = {'mhclovac': []}
frank_scores.update({k: [] for k in IEDB_METHODS})

np.random.seed(0)
random_selection = np.random.randint(0, 1659, N_SAMPLES)

for i, (sequence_name, sequence) in enumerate(read_fasta(bench_fasta)):
    if i not in random_selection:
        continue

    print(i, sequence_name)
    true_epitope = sequence_name.split()[0]
    mhc_allele = mhcnames.normalize_allele_name(sequence_name.split()[1])
    with open('temp.fasta', 'w') as f:
        f.write(f'>{sequence_name}\n{sequence}')
    fasta = os.path.join(os.getcwd(), 'temp.fasta')

    for method in IEDB_METHODS:
        try:
            frank = get_iedb_method_frank(fasta, mhc_allele, true_epitope, method)
            frank_scores[method].append(frank)
        except Exception as e:
            print(method, e)
            frank_scores[method].append(None)
    try:
        frank = get_mhclovac_frank(sequence, mhc_allele, true_epitope)
        frank_scores['mhclovac'].append(frank)
    except Exception as e:
        print('mhclovac', e)
        frank_scores['mhclovac'].append(None)


data = pd.DataFrame.from_dict(frank_scores)
data.to_csv('./results/benchmarking_results.csv', index=False)
