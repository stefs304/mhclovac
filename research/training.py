from mhclovac.training import train_binding_model
from mhclovac.preprocessing import get_features
from mhclovac.config import Config
import pandas as pd
import multiprocessing as mp


TRAINING_SET_SIZE_THRESHOLD = 50
RANDOM_SEED = 0
N_CPU = 8


def train_model(data, mhc_allele):
    x = get_features(peptide_list=data['peptide'], index_id_list=Config.INDEX_ID_LIST)
    y = data['target']
    train_binding_model(features=x, target_values=y, mhc_name=mhc_allele, random_state=RANDOM_SEED)
    return None


train_data = pd.read_csv('../data/ba_el_data.zip')
pool = mp.Pool(processes=N_CPU)

for mhc_key in list(train_data['mhc_allele'].unique()):

    if mhc_key not in ['Patr-A*01:01']:
        continue

    mhc_train_data = train_data[train_data['mhc_allele'] == mhc_key]

    if len(mhc_train_data) < TRAINING_SET_SIZE_THRESHOLD:
        continue

    pool.apply_async(train_model, args=(mhc_train_data, mhc_key))


pool.close()
pool.join()

