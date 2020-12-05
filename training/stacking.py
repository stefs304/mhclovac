import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd
import numpy as np
import pickle
from sklearn.metrics import r2_score
from mhclovac.models import BindingModel


# Setup plot formatting, fonts etc.
rcParams['font.family'] = 'sans-serif'
# rcParams['font.sans-serif'] = ['Tahoma']
rcParams['font.size'] = 12
props = dict(boxstyle='round', facecolor='wheat', alpha=0.6)


def meas_func(x):
    """
    Transform ic50 values into binding score.
    :param x: float; ic50 value
    :return: float
    """
    return (10 - np.log(x)) / 10


def label_func(x):
    """
    Transform qualitative label into binary label.
    :param x: str; label
    :return: int
    """
    postives = [
        'Positive-High',
        'Positive-Intermediate',
        'Positive',
        'Positive-Low'
    ]
    if x in postives:
        return 1
    else:
        return 0


index_id_list = [
    'PRAM900101',
    'FASG760101',
    'ZIMJ680104',
    'CHOP780201',
    'PRAM820103',
    'RACS820112',
    'ROBB760107'
]

train_set_percentage = 0.9

train_set_size_threshold = 50


data = pd.read_csv(f'../data/combined_data.tsv', sep='\t', low_memory=False)
data['meas'] = data['meas'].apply(lambda x: 20000 if x > 20000 else x)  # cap at 20k
data['meas'] = data['meas'].apply(lambda x: 1.0 if x < 1.0 else x)
data['pep_len'] = data['sequence'].apply(len)
# data.drop(index=data[data['pep_len'] < 7].index, axis=0, inplace=True)


for mhc_key in list(data['mhc'].unique()):

    if mhc_key != 'HLA-A*02:01':
        continue

    print(mhc_key)

    tmp_data = data[data['mhc'] == mhc_key]
    tmp_bdata = tmp_data[['sequence', 'meas']].dropna()

    print(f'n samples: {len(tmp_bdata)}')

    if len(tmp_bdata) >= train_set_size_threshold:

        train = tmp_bdata.sample(frac=train_set_percentage)
        test = tmp_bdata.drop(index=train.index)

        X_train = train['sequence']
        y_train = train['meas']

        X_test = test['sequence']
        y_test = test['meas']

        # get the models to evaluate
        model = BindingLovac()
        model.fit(X_train, y_train)

        predictions = model.predict(X_test)
        r2 = round(r2_score(y_test.apply(meas_func), predictions), 3)
        print(f'r2 = {r2}')

        plt.scatter(y_test.apply(meas_func), predictions, s=4)
        plt.savefig('test.png')
