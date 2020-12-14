import subprocess
from mhclovac.sequence import read_fasta, chop_sequence
import mhclovac
import pandas as pd
import mhcnames
import os
from io import StringIO
import seaborn as sns
import matplotlib.pyplot as plt


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

for i, (sequence_name, sequence) in enumerate(read_fasta(bench_fasta)):
    if i > N_SAMPLES:
        break

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


plt.figure(figsize=(10, 6), dpi=100)

sns.set_theme(style="whitegrid")
palette = {k: 'blue' for k in data.columns}
palette['mhclovac'] = 'red'
ax = sns.boxplot(data=data, color='#4c72b0')
ax.artists[0].set_facecolor('red')
ax.set_xticklabels([k for k in data.columns], rotation=30)
ax.set_yscale('log')
ax.set_ylim([0.0005, 1.0])
plt.savefig('figures/mhclovac-frank-benchmark.png', bbox_inches='tight')