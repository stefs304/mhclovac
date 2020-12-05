from .utils import load_model
from .models import BindingModel, LigandModel
from .preprocessing import get_features, transform_ic50_measures, transform_qualitative_measures
from .config import Config
import os
import joblib


def create_model_template(mhc_name):
    model = {
        'mhc_key': mhc_name,
        'binding_model': None,
        'ligand_model': None
    }
    return model


def train_binding_model(peptide_list, ic50_values, mhc_name, verbose=False):
    try:
        model = load_model(mhc_name)
    except:
        model = create_model_template(mhc_name)
    if verbose:
        print(f'Training binding model for {mhc_name}')
        print('Modeling physicochemical properties... This may take a while...')
    X = get_features(peptide_list, Config.INDEX_ID_LIST)
    y = transform_ic50_measures(ic50_values)
    if verbose:
        print(f'n_samples: {X.shape[0]}')
        print(f'n_features: {X.shape[1]}')
    binding_model = BindingModel(verbose=verbose)
    binding_model.fit(X, y)
    model['binding_model'] = binding_model
    dir_path = os.path.dirname(os.path.abspath(__file__))
    filename = os.path.join(dir_path, 'models', f'{mhc_name}.model.gz')
    joblib.dump(model, filename, compress=('gzip', 5))
    if verbose:
        print(f'Model saved; {filename}')
    return None


def train_ligand_model(peptide_list, label_list, mhc_name, verbose=False):
    pass

