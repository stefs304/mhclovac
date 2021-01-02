import pandas as pd
import mhcnames
from mhclovac.sequence import validate_sequence
import os


allele_list = []
with open('/home/stefan/Documents/NetMHCpan_train/allelelist') as f:
    for line in f:
        key = line.split()[0]
        alist = line.split()[1].split(',')
        if len(alist) == 1:
            allele_list.append(alist[0])


def parse_allele_name(allele_name):
    if allele_name not in allele_list:
        return None
    try:
        name = mhcnames.normalize_allele_name(allele_name)
        return name
    except:
        return None


data = []
for file_name in os.listdir('/home/stefan/Documents/NetMHCpan_train'):
    if file_name in ['MHC_pseudo.dat', 'allelelist']:
        continue
    print(file_name)
    df = pd.read_csv(os.path.join('/home/stefan/Documents/NetMHCpan_train', file_name), sep=' ', names=['peptide', 'target', 'mhc_allele'])
    df['mhc_allele'] = df['mhc_allele'].apply(parse_allele_name)
    df.dropna(inplace=True)
    data.append(df)

combined_data = pd.concat(data)

combined_data = combined_data[combined_data['peptide'].apply(lambda x: validate_sequence(x, silent=True))]
combined_data.drop_duplicates(inplace=True)
combined_data.dropna(inplace=True)
print(len(combined_data))

combined_data.to_csv('../data/ba_el_data.zip', index=False)

