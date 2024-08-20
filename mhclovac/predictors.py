
from sklearn.linear_model import LinearRegression
from gaussprot import GaussProt


class BindingPredictor:

    def __init__(self, mhc_allele: str, median_length: int):
        self.model = LinearRegression()
        self.mhc_allele = mhc_allele
        self.median_length = median_length

    def train(self, peptides, alleles, scores):
        pass

    def predict(self, peptides, alleles):
        pass


class PanBindingPredictor:

    def __init__(self, species: str, mhc_class: str, median_length: int):
        self.species = species
        self.mhc_class = mhc_class
        self.median_length = median_length

    def train(self):
        pass

    def predict(self, peptides, alleles):
        pass

