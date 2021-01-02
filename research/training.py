from mhclovac.training import train_binding_model
from mhclovac.preprocessing import get_features
from mhclovac.config import Config
from mhclovac.utils import list_mhc_alleles
import pandas as pd


TRAINING_SET_SIZE_THRESHOLD = 50
RANDOM_SEED = 0

train_data = pd.read_csv('../data/ba_el_data.zip')

for mhc_key in list(train_data['mhc_allele'].unique()):

    if not mhc_key.startswith('Patr'):
        continue

    mhc_train_data = train_data[train_data['mhc_allele'] == mhc_key]

    if len(mhc_train_data) < TRAINING_SET_SIZE_THRESHOLD:
        continue

    print(mhc_key)
    print(f'n_samples = {len(mhc_train_data)}')

    x = get_features(peptide_list=mhc_train_data['peptide'], index_id_list=Config.INDEX_ID_LIST)
    y = mhc_train_data['target']
    train_binding_model(features=x, target_values=y, mhc_name=mhc_key, random_state=RANDOM_SEED)


print(list_mhc_alleles())
