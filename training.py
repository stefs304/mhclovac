import pandas as pd
from mhclovac.training import train_binding_model


TRAINING_SET_SIZE_THRESHOLD = 50


data = pd.read_csv(f'data/combined_data.tsv', sep='\t', low_memory=False)
data['meas'] = data['meas'].apply(lambda x: 20000 if x > 20000 else x)  # cap at 20k


for mhc_key in list(data['mhc'].unique()):
    if mhc_key not in ['H2-Kb']:
        continue
    tmp_data = data[data['mhc'] == mhc_key]
    tmp_bdata = tmp_data[['sequence', 'meas']].dropna().copy()
    if len(tmp_bdata) < TRAINING_SET_SIZE_THRESHOLD:
        continue

    peptide_list = tmp_bdata['sequence']
    ic50_values = tmp_bdata['meas']

    train_binding_model(peptide_list=peptide_list, ic50_values=ic50_values, mhc_name=mhc_key, verbose=True)


