
import os
import pickle
import json
from typing import Union

from sklearn.linear_model import LinearRegression, SGDRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.svm import SVR, LinearSVR
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
pseudoseq_file = os.path.join(os.path.dirname(__file__), 'data', 'pseudo.json')


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
        self.model = MLPRegressor(learning_rate_init=0.01)

    def train(self, peptides, mhcs, values):
        assert len(peptides) == len(mhcs)
        assert len(peptides) == len(values)
        peptide_features = self._generate_peptide_features(peptides)
        mhc_features = self._generate_mhc_features(mhcs)
        features = pd.concat([peptide_features, mhc_features], axis=1)
        self.model.fit(features, values)

    def predict(self, peptides, mhcs: Union[str, list]):
        if isinstance(mhcs, str):
            mhcs = [mhcs for _ in range(len(peptides))]
        else:
            assert len(peptides) == len(mhcs)
        peptide_features = self._generate_peptide_features(peptides)
        mhc_features = self._generate_mhc_features(mhcs)
        features = pd.concat([peptide_features, mhc_features], axis=1)
        return self.model.predict(features)

    def predict_from_pseudosequences(self, peptides, pseudosequences: Union[str, list]):
        raise NotImplementedError()

    def _generate_peptide_features(self, peptides):
        df = pd.DataFrame()
        for schema in schemas:
            core = MhclovacCore(schema, discrete_model_length=self.median_length)
            data = core.generate_models(peptides)
            df = pd.concat([df, pd.DataFrame(data)], axis=1)
        return df

    @staticmethod
    def _generate_mhc_features(mhcs):
        df = pd.DataFrame()
        if not os.path.exists(pseudoseq_file):
            raise FileNotFoundError(f'no pseudosquence data found')
        with open(pseudoseq_file, 'r') as f:
            pseudoseq_dict = json.load(f)
        pseudosequences = [pseudoseq_dict[mhc] for mhc in mhcs]
        for schema in schemas:
            core = MhclovacCore(schema)
            data = core.simply_encode(pseudosequences)
            df = pd.concat([df, pd.DataFrame(data)], axis=1)
        return df


