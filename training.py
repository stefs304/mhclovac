import pandas as pd
from mhclovac.training import train_binding_model


TRAINING_SET_SIZE_THRESHOLD = 50
RANDOM_SEED = 0

data = pd.read_csv(f'data/combined_data.csv')

# print(f'Total number of samples: {data.shape[0]}')


for mhc_key in list(data['mhc'].unique()):

    # if mhc_key not in ['HLA-B*44:02']: continue

    mhc_data = data[data['mhc'] == mhc_key]

    if len(mhc_data) < TRAINING_SET_SIZE_THRESHOLD:
        # skip if less than TRAINING_SET_SIZE_THRESHOLD peptides
        continue

    peptide_list = mhc_data['peptide']
    ic50_values = mhc_data['ic50']

    train_binding_model(
        peptide_list=peptide_list,
        ic50_values=ic50_values,
        mhc_name=mhc_key,
        verbose=2,
        random_state=RANDOM_SEED
    )

