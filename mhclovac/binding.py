
import os
import pickle
import json
from sklearn.linear_model import LinearRegression
from sklearn.svm import SVR
import pandas as pd
from mhclovac.core import MhclovacCore
from mhclovac.const import (
    hydrogen_bonds_params,
    wimley_white_params,
    polarizability_params,
    isoelectric_params,
    class_1_pep_len,
    class_2_pep_len
)


models_file = os.path.join(os.path.dirname(__file__), 'data', 'trained_models.pickle')
pseudoseq_file = os.path.join(os.path.dirname(__file__), 'data', 'pseudoseq.txt')

schemas = [
    hydrogen_bonds_params,
    wimley_white_params,
    polarizability_params,
]

class MhcSpecificBindingPredictor:

    def __init__(self, mhc: str, class_: int):
        self.model = SVR()
        self.mhc = mhc
        self.median_length = class_2_pep_len if class_ == 2 else class_1_pep_len

    def train(self, peptides, values):
        assert len(peptides) == len(values)
        features = self._get_features(peptides)
        self.model.fit(features, values)

    def predict(self, peptides):
        features = self._get_features(peptides)
        return self.model.predict(features)

    def _get_features(self, peptides):
        df = pd.DataFrame()
        for schema in schemas:
            core = MhclovacCore(schema, discrete_model_length=self.median_length)
            data = core.generate_models(peptides)
            df = pd.concat([df, pd.DataFrame(data)], axis=1)
        return df

    @classmethod
    def load(cls, mhc):
        with open(models_file, 'rb') as f:
            models = pickle.load(f)
            try:
                return models[mhc]
            except KeyError:
                raise KeyError(f'Unsupported mhc: {mhc}')




class PanSpecificBindingPredictor:

    def __init__(self, class_: int):
        self.median_length = class_2_pep_len if class_ == 2 else class_1_pep_len
        self.model = SVR()

    def train(self, peptides, alleles, values):
        raise NotImplementedError()

    def predict(self, peptides, allele):
        raise NotImplementedError()

    def load(self):
        pass

    def _generate_peptide_features(self, peptides):
        pass

    def _generate_allele_features(self, alleles):
        pass



