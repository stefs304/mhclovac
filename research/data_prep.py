import pandas as pd
import mhcnames
from mhclovac.sequence import validate_sequence


def parse_allele_name(allele_name):
    try:
        name = mhcnames.normalize_allele_name(allele_name)
        return name
    except:
        return None


binding_data_file = '~/Documents/bdata.20130222.mhci.txt'
ligand_data_file = '~/Documents/mhc_ligand_full.csv'

binding_data = pd.read_csv(
    binding_data_file,
    sep='\t',
    skiprows=1,
    usecols=[1, 3, 5],
    names=['mhc_allele', 'peptide', 'affinity']
)
binding_data.dropna(inplace=True)
print(f'BA data, n_samples={len(binding_data)}')


# 10 - Linear peptide
# 11 - Description - peptide
# 80 - Assay Group
# 81 - Units
# 85 - Quantitative measurement
# 83 - Qualitative Measure
# 95 - Allele Name
# 98 - MHC allele class - I
ligand_data = pd.read_csv(
    ligand_data_file,
    sep=',',
    skiprows=2,
    nrows=None,
    usecols=[10, 11, 80, 81, 85, 95, 98],
    names=['peptide_type', 'peptide', 'assay_group', 'units', 'affinity', 'mhc_allele', 'mhc_class']
)

ligand_data = ligand_data[ligand_data['mhc_class'] == 'I']

assay_group_list = [
    # 'dissociation constant KD (~EC50)',
    'half maximal inhibitory concentration (IC50)',
    'dissociation constant KD (~IC50)',
    # 'half maximal effective concentration (EC50)',
    'dissociation constant KD'
]

ligand_data = ligand_data[ligand_data['assay_group'].isin(assay_group_list)]
ligand_data = ligand_data[['peptide', 'affinity', 'mhc_allele']]
ligand_data.dropna(inplace=True)

ligand_data['mhc_allele'] = ligand_data['mhc_allele'].apply(parse_allele_name)
ligand_data.dropna(inplace=True)
print(f'EL data, n_samples={len(ligand_data)}')

combined_data = pd.concat([binding_data, ligand_data])
combined_data = combined_data[combined_data['peptide'].apply(lambda x: validate_sequence(x, silent=True))]
combined_data.drop_duplicates(inplace=True)
combined_data.dropna(inplace=True)
print(f'Combined data, n_samples={len(combined_data)}')

combined_data.to_csv('../data/combined_data.zip', index=False)
