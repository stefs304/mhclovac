import numpy as np
import pandas as pd
from mhclovac.sequence import model_distribution, encode_sequence
from .utils import load_index_data
import math


def standardize_index(index: dict) -> dict:
    """
    Standardize values of given index to range [-1 1].
    """
    values = np.array([index[k] for k in index])
    mu = values.mean()
    std = values.std()
    for k, v in index.items():
        index[k] = (v - mu) / std
    return index


def normalize_index_data(index: dict) -> dict:
    normalized_index = dict(index)
    values = np.array([index[k] for k in index])
    min_ = values.min()
    max_ = values.max()
    for k, v in index.items():
        normalized_index[k] = 2 * (v - min_) / (max_ - min_) - 1
    return normalized_index


def sequence_to_features(sequence: str, index_list: list, sigma: float = 0.4, n_discrete_points: int = 10) -> list:
    """
    Convert sequence string to list of features.
    """
    sequence_features = []
    for index in index_list:
        features = model_distribution(sequence, index, sigma=sigma, n_discrete_points=n_discrete_points)
        sequence_features.extend(features)
    return sequence_features


def anchors_to_features(sequence: str, index_data_list: list):
    if len(sequence) < 4:
        raise ValueError(f'peptide must be longer than 4')
    anchor_features = []
    for index_data in index_data_list:
        features = encode_sequence(sequence=sequence, encoding_scheme=index_data)
        anchor1 = features[:4]
        anchor2 = features[-4:]
        anchor_features.extend(anchor1)
        anchor_features.extend(anchor2)
    return anchor_features


def get_sequence_features(peptide_list, index_id_list):
    index_data = []
    index_data_full = load_index_data()
    peptide_df = pd.DataFrame()
    peptide_df['peptide'] = peptide_list
    for index_id in index_id_list:
        index_data.append(index_data_full[index_id]['index_data'])
    features = peptide_df['peptide'].apply(lambda x: sequence_to_features(x, index_data))
    return pd.DataFrame(features.tolist())


def get_anchor_features(peptide_list, index_id_list):
    index_data = []
    index_data_full = load_index_data()
    peptide_df = pd.DataFrame()
    peptide_df['peptide'] = peptide_list
    for index_id in index_id_list:
        index_data.append(index_data_full[index_id]['index_data'])
    features = peptide_df['peptide'].apply(lambda x: anchors_to_features(x, index_data))
    return pd.DataFrame(features.tolist())


def ic50_to_score(ic50_list):
    y = pd.DataFrame()
    y['ic50'] = ic50_list
    y['ic50'] = y['ic50'].apply(lambda x: (10 - np.log(x)) / 10)
    return y['ic50']
