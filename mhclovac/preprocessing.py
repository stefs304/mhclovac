import numpy as np
import pandas as pd
from mhclovac.sequence import model_distribution
from .utils import load_index_data
from .config import Config
import math
from pandarallel import pandarallel


def normalize_index_data(index: dict) -> dict:
    normalized_index = dict(index)
    values = np.array([index[k] for k in index])
    min_ = values.min()
    max_ = values.max()
    for k, v in index.items():
        normalized_index[k] = 2 * (v - min_) / (max_ - min_) - 1
    return normalized_index


def standardize_index_data(index: dict) -> dict:
    standardized_index = dict(index)
    values = np.array([index[k] for k in index])
    mu = np.mean(values)
    std = np.std(values)
    for k, v in index.items():
        standardized_index[k] = (v - mu) / std
    return standardized_index


def sequence_to_features(sequence: str, index_list: list) -> list:
    """
    Convert sequence string to list of features.
    """
    sequence_features = []
    for index in index_list:
        features = model_distribution(
            sequence=sequence,
            encoding_scheme=index,
            sigma=Config.SIGMA,
            overlap_distance=Config.OVERLAP_DISTANCE,
            n_discrete_points=Config.N_DISCRETE_POINTS
        )
        sequence_features.extend(features)
    return sequence_features


def get_features(peptide_list, index_id_list, n_cpu=1):
    if n_cpu != 1:
        pandarallel.initialize(nb_workers=n_cpu, verbose=0)
    index_data = load_index_data(index_id_list=index_id_list)
    peptide_df = pd.DataFrame()
    peptide_df['peptide'] = peptide_list
    if n_cpu == 1:
        features = peptide_df['peptide'].apply(lambda x: sequence_to_features(x, index_data))
    else:
        features = peptide_df['peptide'].parallel_apply(lambda x: sequence_to_features(x, index_data))
    return pd.DataFrame(features.tolist())


def transform_ic50_values(ic50_values):
    data = pd.DataFrame()
    data['values'] = ic50_values
    return data['values'].apply(lambda x: 1 - math.log(x, 50000) if x > 0 else 0)
