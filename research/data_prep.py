import pandas as pd
import mhcnames
from mhclovac.sequence import validate_sequence


def parse_allele_name(allele_name):
    try:
        name = mhcnames.normalize_allele_name(allele_name)
        return name
    except:
        return None


nrows = None

ligand_data = pd.read_csv('~/Documents/mhc_ligand_full.csv', skiprows=2, low_memory=True, nrows=nrows, usecols=[11, 85, 95, 98], names=['peptide', 'ic50', 'mhc', 'class'])
binding_data = pd.read_csv('~/Documents/bdata.20130222.mhci.txt', sep='\t', skiprows=1, low_memory=False, nrows=nrows, usecols=[1, 3, 5], names=['mhc', 'peptide', 'ic50'])


print(f'original ligand_data n_samples: {ligand_data.shape[0]}')
print(f'original binding_data n_samples: {binding_data.shape[0]}')


ligand_data = ligand_data[ligand_data['class'] == 'I']
ligand_data.drop(columns=['class'], inplace=True)
ligand_data['mhc'] = ligand_data['mhc'].apply(parse_allele_name)
ligand_data.dropna(inplace=True)

binding_data['mhc'] = binding_data['mhc'].apply(parse_allele_name)
binding_data.dropna(inplace=True)

print(f'ligand_data n_samples: {ligand_data.shape[0]}')
print(f'binding_data n_samples: {binding_data.shape[0]}')

combined_data = pd.concat([ligand_data, binding_data])
combined_data['valid_peptide'] = combined_data['peptide'].apply(lambda x: validate_sequence(x, silent=True))
combined_data = combined_data[combined_data['valid_peptide'] == True]
combined_data.drop(columns=['valid_peptide'], inplace=True)
combined_data.drop(index=combined_data[combined_data['ic50'] == 0].index, inplace=True)
# combined_data['ic50'] = combined_data['ic50'].apply(lambda x: 0.1 if x == 0 else x)
combined_data.drop_duplicates(inplace=True)
combined_data.dropna(inplace=True)
print(f'combined_data n_samples: {combined_data.shape[0]}')

combined_data.to_csv('../data/combined_data.csv', index=False)
