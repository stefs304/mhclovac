from .utils import load_model
from .models import BindingModel
from .preprocessing import get_features, transform_ic50_values
from .config import Config
import os
import joblib
import time


def create_model_template(mhc_name):
    model = {
        'mhc_key': mhc_name,
        'binding_model': None
    }
    return model


def train_binding_model(peptide_list, ic50_values, mhc_name, verbose=False, random_state=None):
    try:
        model = load_model(mhc_name)
    except:
        model = create_model_template(mhc_name)
    if verbose:
        print(f'Training binding model for {mhc_name}')
        print(f'n_samples: {len(peptide_list)}')
        print('Modeling physicochemical properties... This may take a while...')
    t0 = time.time()
    X = get_features(peptide_list, Config.INDEX_ID_LIST)
    # y = transform_ic50_values(ic50_values)
    y = ic50_values
    t_elapsed = round(time.time() - t0, 1)
    if verbose:
        print(f'features: {X.shape[0]} x {X.shape[1]} in {t_elapsed} seconds.')
    binding_model = BindingModel(verbose=verbose, random_state=random_state)
    binding_model.fit(X, y)
    model['binding_model'] = binding_model
    dir_path = os.path.dirname(os.path.abspath(__file__))
    filename = os.path.join(dir_path, 'models', f'{mhc_name}.model.gz')
    joblib.dump(model, filename, compress=('gzip', 5))
    if verbose:
        print(f'Model saved; {filename}')
    return None
