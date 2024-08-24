
from sklearn.linear_model import LinearRegression
from sklearn.svm import SVR
import pandas as pd
from mhclovac.core import MhclovacCore
from mhclovac.schemas import (
    hydrogen_bonds_params,
    wimley_white_params,
    polarizability_params,
    isoelectric_params
)


schemas = [
    hydrogen_bonds_params,
    wimley_white_params,
    polarizability_params,
]

class AlleleSpecificBindingPredictor:

    def __init__(self, mhc_allele: str, median_length: int):
        self.model = SVR()
        self.mhc_allele = mhc_allele
        self.median_length = median_length

    def train(self, peptides, values):
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



class PanBindingPredictor:

    def __init__(self, species: str, mhc_class: str, median_length: int):
        self.species = species
        self.mhc_class = mhc_class
        self.median_length = median_length

    def train(self):
        pass

    def predict(self, peptides, alleles):
        pass

