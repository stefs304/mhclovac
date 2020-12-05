import numpy as np
import os
import joblib


def list_mhc_alleles():
    """
    Returns the list of available MHC alleles.

    :return: list
    """
    mhc_alleles = []
    for fname in os.listdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'models')):
        mhc_name = fname.split('.model.gz')[0]
        mhc_alleles.append(mhc_name)

    return mhc_alleles


def load_model(mhc):
    """
    Returns trained models for given MHC allele.

    :param mhc: str, MHC allele
    :return: mhclovac.BindingPredictor, mhclovac.EpitopePredictor
    """
    if mhc not in list_mhc_alleles():
        msg = f'"{mhc}" not available.'
        raise ValueError(msg)
    fpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'models', f'{mhc}.model.gz')
    data = joblib.load(fpath)
    bmodel, emodel = data['binding'], data['epitope']

    return bmodel, emodel


def load_index():
    """
    Loads and returns index data.

    :return: list, index data
    """
    index_list = ['PRAM900101', 'FASG760101', 'ZIMJ680104', 'CHOP780201', 'PRAM820103', 'RACS820112', 'ROBB760107']
    index_data = []
    for index_id in index_list:
        fpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'index', f'{index_id}.index.gz')
        data = joblib.load(fpath)
        index_data.append(data['index_data'])

    return index_data


def load_index_data():
    fpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'index', 'index_data.gz')
    data = joblib.load(fpath)
    return data


def index_correlation(index1, index2):
    """
    Find correlation coefficient for two physicochemical indexes.
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
