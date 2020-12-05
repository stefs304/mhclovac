from .utils import load_model
from .models import BindingModel, LigandModel
from .preprocessing import get_features
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
    binding_model = BindingModel(verbose=verbose)
    binding_model.fit(peptide_list=peptide_list, ic50_values=ic50_values)
    model['binding_model'] = binding_model
    dir_path = os.path.dirname(os.path.abspath(__file__))
    filename = os.path.join(dir_path, 'models', f'{mhc_name}.model.gz')
    joblib.dump(model, filename, compress=('gzip', 5))
    return


def train_ligand_model(peptide_list, label_list, mhc_name, verbose=False):
    pass

