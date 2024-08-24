
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


schemas = [
    hydrogen_bonds_params,
    wimley_white_params,
    polarizability_params,
]

class AlleleSpecificBindingPredictor:

    def __init__(self, mhc_allele: str, class_: int):
        self.model = SVR()
        self.mhc_allele = mhc_allele
        self.median_length = class_2_pep_len if class_ == 2 else class_1_pep_len

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



class PanSpecificBindingPredictor:

    def __init__(self):
        raise NotImplementedError()

    def train(self, peptides, alleles, values):
        raise NotImplementedError()

    def predict(self, peptides, allele):
        raise NotImplementedError()


