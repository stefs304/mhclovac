import pandas as pd
from mhclovac.sequence import validate_sequence


species_map = {
    'human': ['human', 'human (Homo sapiens)', 'Homo sapiens'],
    'mouse': ['mouse (Mus musculus)', 'Mus musculus', 'mouse'],
    'rat': ['rat (Rattus norvegicus)', 'rat', 'Rattus norvegicus'],
    'macaque': ['rhesus macaque (Macaca mulatta)', 'macaque'],
    'gorilla': ['gorilla (Gorilla gorilla)', 'gorilla'],
    'pig': ['pig (Sus scrofa)', 'pig', 'Sus scrofa'],
    'cow': ['cattle (Bos taurus)', 'cow'],
    'horse': ['horse', 'horse (Equus caballus)', 'Equus caballus'],
    'chimpanzee': ['chimpanzee', 'chimpanzee (Pan troglodytes)', 'Pan troglodytes'],
    'chicken': ['chicken (Gallus gallus)', 'Gallus gallus'],
    'dog': ['dog (Canis lupus familiaris)', 'Canis lupus familiaris'],
    'cat': ['cat (Felis catus)']
}


def get_species(label):
    for key in species_map:
        if label in species_map[key]:
            return key
    None


mhc_map = {
    'H-2-Db': 'H2-Db',
    'H-2-Dd': 'H2-Dd',
    'H-2-Kb': 'H2-Kb',
    'H-2-Kd': 'H2-Kd',
    'H-2-Kk': 'H2-Kk',
    'H-2-Ld': 'H2-Ld',
    'H-2-Lq': 'H2-Lq',
    'RT1A': 'RT1-A',
    'SLA-1*0401': 'SLA-1*04:01',
    'SLA-1*0701': 'SLA-1*07:01',
    'SLA-2*0401': 'SLA-2*04:01',
    'SLA-3*0401': 'SLA-3*04:01',
    'ELA-A1': 'ELA-A1 class I',
    'Patr-A*0101': 'Patr-A*01:01',
    'Patr-A*0301': 'Patr-A*03:01',
    'Patr-A*0401': 'Patr-A*04:01',
    'Patr-A*0602': 'Patr-A*06:02',
    'Patr-A*0701': 'Patr-A*07:01',
    'Patr-A*0901': 'Patr-A*09:01',
    'Patr-B*0101': 'Patr-B*01:01',
    'Patr-B*0901': 'Patr-B*09:01',
    'Patr-B*1301': 'Patr-B*13:01',
    'Patr-B*1701': 'Patr-B*17:01',
    'Patr-B*2401': 'Patr-B*24:01',

}

bdata = pd.read_csv('/home/stefan/Documents/bdata.20130222.mhci.txt', sep='\t', low_memory=False)
edata = pd.read_csv('/home/stefan/Documents/edata_sstojanovic.tsv', sep='\t')

for species in species_map:

    print(f'{species}')

    tmp_bdata = bdata[bdata['species'].isin(species_map[species])]
    tmp_edata = edata[edata['species'].isin(species_map[species])]

    bdata_mhc_list = list(tmp_bdata['mhc'].unique())
    edata_mhc_list = list(tmp_edata['mhc'].unique())

    non_match_mhc_list = []

    print(f'\tin bdata, not in edata')
    for mhc_key in bdata_mhc_list:
        if mhc_key not in edata_mhc_list:
            print(f'\t{mhc_key}')

    print(f'\tin edata, not in bdata')
    for mhc_key in edata_mhc_list:
        if mhc_key not in bdata_mhc_list:
            print(f'\t{mhc_key}')


bdata['mhc'] = bdata['mhc'].apply(lambda x: mhc_map[x] if x in mhc_map else x)

data = edata.append(bdata)

print(data.shape)

data['meas'] = data['meas'].apply(lambda x: 0.1 if x == 0 else x)
data['meas'] = data['meas'].apply(lambda x: 20000 if x > 20000 else x)

data['valid_seq'] = data['sequence'].apply(lambda x: validate_sequence(x))
data = data[data['valid_seq'] == True]
data.drop(['valid_seq'], axis=1, inplace=True)

print(data.shape)

data['species'] = data['species'].apply(get_species)
data = data[~data['species'].isna()]

print(data.shape)

print(data.shape)
print(data['species'].unique())

data.to_csv('../data/combined_data.tsv', sep='\t')
