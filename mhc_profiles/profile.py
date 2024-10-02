
import os
import pickle
import numpy as np
from dataclasses import dataclass
from .constants import peptide_length_map
from .core import PhysicochemicalProfiler


file_dir = os.path.abspath(os.path.dirname(__file__))
static_dir = os.path.join(file_dir, "static")
data_file = os.path.join(static_dir, "data.pickle")


class MHCProfile:

    def __init__(self, mhc, species, schema_name, schema_data, mhc_class):
        self.mhc = mhc
        self.species = species
        self.mhc_class = mhc_class
        self.schema_name = schema_name
        self.schema_data = schema_data
        self.median_peptide_length = peptide_length_map[mhc_class]
        self.components = [PositionalComponent() for _ in range(self.median_peptide_length)]
        self.pcp = PhysicochemicalProfiler(schema=self.schema_data, discrete_model_length=self.median_peptide_length)

    def train(self, peptides, reference_peptides):
        profiles = self.pcp.generate_models(peptides)
        ref_profiles = self.pcp.generate_models(reference_peptides)
        for i in range(self.median_peptide_length):
            self.components[i].compute_stats(profiles[:, i], ref_profiles[:, i])

    def predict(self, peptides):
        profiles = self.pcp.generate_models(peptides)
        pass

    @property
    def profile_stats(self):
        return None


class PositionalComponent:

    def __init__(self):
        self._mu = None
        self._std = None
        self._ref_mu = None
        self._ref_std = None

    def compute_stats(self, x, x_ref):
        self._mu = np.mean(x)
        self._std = np.std(x)
        self._ref_mu = np.mean(x_ref)
        self._ref_std = np.std(x_ref)

    def evaluate(self, x):
        pass


@dataclass
class MHCProfileData:
    mhc: str
    profiles: dict[str, MHCProfile]


class MHCI:

    def __init__(self, profile_data):
        self._profile_data = profile_data

    @staticmethod
    def load(mhc: str):
        with open(data_file, "rb") as f:
            data = pickle.load(f)
        if mhc in data:
            return MHCI(data[mhc])
        raise ValueError(f'mhc not found in data')

    def save(self, overwrite=False):
        if not os.path.exists(data_file):
            with open(data_file, "wb") as f:
                pickle.dump({}, f)
        with open(data_file, "rb") as f:
            data = pickle.load(f)
        if self.mhc in data and not overwrite:
            raise ValueError("MHC already exists.")
        data[self.mhc] = self._profile_data
        with open(data_file, "wb") as f:
            pickle.dump(data, f)

    @property
    def profiles(self) -> dict:
        return self._profile_data.profiles

    @property
    def mhc(self) -> str:
        return self._profile_data.mhc

    def predict(self):
        pass


