import pandas as pd
from mhclovac.training import train_binding_model
import multiprocessing as mp


TRAINING_SET_SIZE_THRESHOLD = 50
RANDOM_SEED = 0
N_PROC = 2

data = pd.read_csv(f'data/train_data.zip')

pool = mp.Pool(processes=N_PROC)

for mhc_key in list(data['mhc'].unique()):

    if mhc_key not in ['HLA-B*44:02']: continue

    mhc_data = data[data['mhc'] == mhc_key]

    if len(mhc_data) < TRAINING_SET_SIZE_THRESHOLD:
        # skip if less than TRAINING_SET_SIZE_THRESHOLD peptides
        continue

    peptide_list = mhc_data['peptide']
    target_values = mhc_data['target']

    pool.apply_async(train_binding_model, (peptide_list, target_values, mhc_key, False, RANDOM_SEED, 1))
    # args: peptide_list, target_values, mhc_name, verbose, random_seed, n_jobs
pool.close()
pool.join()
