import numpy as np
import pandas as pd
from mhclovac.sequence import model_distribution
from .utils import load_index_data
import math


def normalize_index_data(index: dict) -> dict:
    normalized_index = dict(index)
    values = np.array([index[k] for k in index])
    min_ = values.min()
    max_ = values.max()
    for k, v in index.items():
        normalized_index[k] = 2 * (v - min_) / (max_ - min_) - 1
    return normalized_index


def sequence_to_features(sequence: str, index_list: list, n_discrete_points: int = 10) -> list:
    """
    Convert sequence string to list of features.
    """
    sequence_features = []
    for index in index_list:
        features = model_distribution(sequence, index, n_discrete_points=n_discrete_points)
        sequence_features.extend(features)
    return sequence_features


def get_features(peptide_list, index_id_list):
    index_data = load_index_data(index_id_list=index_id_list)
    peptide_df = pd.DataFrame()
    peptide_df['peptide'] = peptide_list
    features = peptide_df['peptide'].apply(lambda x: sequence_to_features(x, index_data))
    return pd.DataFrame(features.tolist())


def get_label(measure: str) -> int:
    """
    Transform qualitative measure.
    """
    positive_measures = [
        'Positive-High',
        'Positive-Intermediate',
        'Positive',
        'Positive-Low'
    ]
    if measure in positive_measures:
        return 1
    return 0


def transform_ic50_measures(ic50_values):
    data = pd.DataFrame()
    data['values'] = ic50_values
    return data['values'].apply(lambda x: 10 - np.log(x) / 10)


def transform_qualitative_measures(value_list):
    data = pd.DataFrame()
    data['values'] = value_list
    return data['values'].apply(get_label)
