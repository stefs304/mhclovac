from .utils import save_trained_model
from .models import BindingModel
import time


def train_binding_model(features, target_values, mhc_name, verbose=False, random_state=None, n_jobs=-1):
    if verbose:
        print(f'Training binding model for {mhc_name}')
        print(f'n_samples: {len(features)}')
        print('Modeling physicochemical properties... This may take a while...')
    t0 = time.time()
    t_elapsed = round(time.time() - t0, 1)
    if verbose:
        print(f'features: {features.shape[0]} x {features.shape[1]} in {t_elapsed} seconds.')
    binding_model = BindingModel(verbose=verbose, random_state=random_state, n_jobs=n_jobs)
    binding_model.fit(features, target_values)
    save_trained_model(mhc_allele=mhc_name, trained_model=binding_model)
    if verbose:
        print(f'Trained {mhc_name} model saved.')
    return None
