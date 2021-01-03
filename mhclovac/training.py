from .models import BindingModel
import os
import joblib


def train_binding_model(features, target_values, mhc_name, verbose=False, random_state=None, n_jobs=-1):
    """
    Train binding model and save file.
    """
    model = BindingModel(verbose=verbose, random_state=random_state, n_jobs=n_jobs)
    model.fit(features, target_values)
    dir_path = os.path.dirname(os.path.abspath(__file__))
    filename = os.path.join(dir_path, 'models', f'{mhc_name}.model.gz')
    joblib.dump(model, filename, compress=('gzip', 5))
    if verbose:
        print(f'Model saved; {filename}')
    return None
