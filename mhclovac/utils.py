import numpy as np
import os
import joblib


def list_mhc_alleles():
    """
    Returns the list of available MHC alleles.
    """
    mhc_alleles = []
    for fname in os.listdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'models')):
        mhc_name = fname.split('.model.gz')[0]
        mhc_alleles.append(mhc_name)
    return mhc_alleles


def load_model(mhc: str):
    """
    Load trained models for a given MHC allele.
    """
    if mhc not in list_mhc_alleles():
        msg = f'"{mhc}" not supported.'
        raise ValueError(msg)
    model_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'models', f'{mhc}.model.gz')
    model = joblib.load(model_path)
    return model


def load_index_data(index_id_list: list = None) -> list or dict:
    """
    Load index data for a given list of index IDs.
    """
    file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'index', 'index_data.gz')
    data = joblib.load(file_path)
    if index_id_list:
        index_data = []
        for index_id in index_id_list:
            index_data.append(data[index_id]['standardized_index_data'])
        return index_data
    return data


def get_index_correlation(index1, index2):
    """
    Get correlation coefficient for two physicochemical indexes.
    """
    keys = list(index1.keys())
    array_1 = []
    array_2 = []
    for key in keys:
        array_1.append(index1[key])
        array_2.append(index2[key])
    return np.corrcoef(array_1, array_2)[0][1]


def pdf(x, sigma):
    """
    Calculates normal probability density function at x data points. Assumes
    mean of 0 and std provided by sigma parameter.
    :param x: data points
    :param sigma: std
    :return: np.array
    """
    y = np.exp(-x**2/(2*sigma))/(sigma*np.sqrt(2*np.pi))
    y = (y - y.min())/(y.max()-y.min())
    return y
