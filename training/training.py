import os
import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from mhclovac.preprocessing import sequence_to_features
from mhclovac.models import BindingPredictor, EpitopePredictor
from matplotlib import rcParams
from sklearn.metrics import f1_score, r2_score
import joblib
import sys
import json


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


index_id_list = ['PRAM900101', 'FASG760101', 'ZIMJ680104', 'CHOP780201', 'PRAM820103', 'RACS820112', 'ROBB760107']
# PRAM900101 Hydrophobicity (Prabhakaran, 1990) 34.336333333333314
# FASG760101 Molecular weight (Fasman, 1976) 32.02166666666666
# ZIMJ680104 Isoelectric point (Zimmerman et al., 1968) 30.74800000000001
# CHOP780201 Normalized frequency of alpha-helix (Chou-Fasman, 1978b) 27.77033333333333
# PRAM820103 Correlation coefficient in regression analysis (Prabhakaran-Ponnuswamy, 1982) 26.523333333333337
# RACS820112 Average relative fractional occurrence in ER(i-1) (Rackovsky-Scheraga, 1982) 22.355000000000004
# ROBB760107 Information measure for extended without H-bond (Robson-Suzuki, 1976) 19.51566666666668


train_set_percentage = 0.9

train_set_size_threshold = 50

top_score_threshold = 20

score_percentile_threshold = 90

correlation_threshold = 0.7

skip_list = [
    'HLA class I',
    'B19 class I',
    'H2-Kb Y22F, M23I, E24S, D30N mutant',
    'RT1-A',
]

with open('../data/aaindex_data.pickle', 'rb') as f:
    aaindex_data = pickle.load(f)

index_data = []
for id_ in index_id_list:
    index_data.append(aaindex_data[id_]['index_data'])


# Load trained models, use to skip alleles that were already processed
trained_models = []
for fname in os.listdir('./trained_models'):
    mhc = fname.split('.model.gz')[0]
    trained_models.append('.'.join(mhc))



data = pd.read_csv(f'../data/combined_data.tsv', sep='\t')
data['meas'] = data['meas'].apply(lambda x: 20000 if x > 20000 else x)  # cap at 20k


for mhc_key in list(data['mhc'].unique()):

    # if mhc_key in trained_models:
    #     print(f'skipping {mhc_key}')
    #     continue

    if mhc_key in skip_list:
        print(f'skipping {mhc_key}')
        continue

    # if mhc_key != 'H2-Kb':
    #     continue

    print(mhc_key)

    models = {
        'mhc_key': mhc_key,
        'binding': None,
        'epitope': None,
        'r2_score': None,
        'f1_score': None
    }
    config = {
        'mhc_key': mhc_key,
        'index_id_list': None
    }

    tmp_data = data[data['mhc'] == mhc_key]

    print(f'data shape: {tmp_data.shape}')

    index_data = []
    for index_id in index_id_list:
        index_data.append(aaindex_data[index_id]['index_data'])

    fig, ax = plt.subplots(2,2, figsize=(10,10))
    fig.suptitle(f'{mhc_key} training stats', fontsize=14)

    tmp_bdata = tmp_data[['sequence', 'meas']].dropna().copy()
    if len(tmp_bdata) >= train_set_size_threshold:

        print('Training binding prediction')

        train = tmp_bdata.sample(frac=train_set_percentage)
        test = tmp_bdata.drop(index=train.index)

        # train = tmp_bdata.sample(frac=1.0)
        # test = train.copy()

        X_train = train['sequence'].apply(lambda x: sequence_to_features(x, index_data))
        X_train = pd.DataFrame(X_train.tolist(), index=train.index)
        y_train = train['meas'].apply(meas_func)

        ax[0][0].hist(y_train, bins=50)
        ax[0][0].grid()
        ax[0][0].set_xlabel('hist of transformed ic50 values')
        ax[0][0].text(0.05, 0.95, f'n_samples = {X_train.shape[0]}', transform=ax[0][0].transAxes, fontsize=12,
                      verticalalignment='top', bbox=props)

        print('model fit')
        # ----- setup model
        try:
            model = BindingPredictor(random_state=None, standardize_data=True)
            print(f'estimators: {model.estimators}')
            model.fit(X_train, y_train)

            print('evaluating model')
            X_test = test['sequence'].apply(lambda x: sequence_to_features(x, index_data))
            X_test = pd.DataFrame(X_test.tolist(), index=test.index)
            y_test = test['meas'].apply(meas_func)

            predictions = model.predict(X_test)
            r2 = round(r2_score(y_test, predictions), 3)
            print(f'r2 = {r2}')

            ax[0][1].scatter(y_test, predictions, s=4)
            ax[0][1].set_xlabel('x = (10 - log(ic50)) / 10 ; (higher the better)')
            ax[0][1].set_ylabel('predicted score')
            # ax[0][1].set_xlim([y_test.min(), y_test.max()])
            # ax[0][1].set_ylim([y_test.min(), y_test.max()])
            ax[0][1].grid()
            ax[0][1].text(0.05, 0.95, f'r2 = {r2}', transform=ax[0][1].transAxes, fontsize=12,
                          verticalalignment='top', bbox=props)

            # save model
            models['binding'] = model
            models['r2_score'] = r2
        except Exception as e:
            msg = f'Error at {mhc_key} binding: {e}\n'
            sys.stderr.write(msg)
            pass
    else:
        ax[0][0].text(0.05, 0.95, f'No data', transform=ax[0][0].transAxes, fontsize=12,
                      verticalalignment='top', bbox=props)
        ax[0][1].text(0.05, 0.95, f'No data', transform=ax[0][1].transAxes, fontsize=12,
                      verticalalignment='top', bbox=props)

    # Now the same for epitopes
    tmp_edata = tmp_data[['sequence', 'label']].dropna().copy()

    if len(tmp_edata) >= train_set_size_threshold:

        print('Training epitope prediction')

        train = tmp_edata.sample(frac=train_set_percentage)
        test = tmp_edata.drop(index=train.index)

        # train = tmp_edata.sample(frac=1.0)
        # test = train.copy()

        X_train = train['sequence'].apply(lambda x: sequence_to_features(x, index_data))
        X_train = pd.DataFrame(X_train.tolist(), index=train.index)
        y_train = train['label'].apply(label_func)

        y_train.value_counts().plot(kind='bar', ax=ax[1][0])
        ax[1][0].grid()
        ax[1][0].set_xlabel('label counts')
        ax[1][0].text(0.05, 0.95, f'n_samples = {X_train.shape[0]}', transform=ax[1][0].transAxes, fontsize=12,
                      verticalalignment='top', bbox=props)

        try:
            print('model fit')
            model = EpitopePredictor(random_state=None, standardize_data=True)
            print(f'estimators: {model.estimators}')
            model.fit(X_train, y_train)

            print('evaluating model')
            X_test = test['sequence'].apply(lambda x: sequence_to_features(x, index_data))
            X_test = pd.DataFrame(X_test.tolist(), index=test.index)
            y_test = test['label'].apply(label_func)

            predictions = model.predict_proba(X_test)
            preds_binary = model.predict(X_test)
            f1 = round(f1_score(y_test, preds_binary), 3)
            print(f'f1 = {f1}')

            by_cat_preds = pd.DataFrame()
            by_cat_preds['label'] = y_test
            by_cat_preds['score'] = predictions

            g = sns.stripplot(x="label", y="score", data=by_cat_preds, ax=ax[1][1], s=4)
            ax[1][1].set_xlabel('labels')
            ax[1][1].set_ylabel('epitope score')
            ax[1][1].grid()
            ax[1][1].text(0.05, 0.95, f'f1 = {f1}', transform=ax[1][1].transAxes, fontsize=12,
                    verticalalignment='top', bbox=props)

            # save model
            models['epitope'] = model
            models['f1_score'] = f1
        except Exception as e:
            msg = f'Error at {mhc_key} epitope: {e}\n'
            sys.stderr.write(msg)
            pass
    else:
        ax[1][0].text(0.05, 0.95, f'No data', transform=ax[1][0].transAxes, fontsize=12,
                      verticalalignment='top', bbox=props)
        ax[1][1].text(0.05, 0.95, f'No data', transform=ax[1][1].transAxes, fontsize=12,
                      verticalalignment='top', bbox=props)

    if '/' in mhc_key:
        mhc_key = mhc_key.replace('/', '-')

    if any([models['binding'], models['epitope']]):

        plt.savefig(f'./results/{mhc_key}.png', bbox_inches='tight')
        joblib.dump(models, f'trained_models/{mhc_key}.model.gz', compress=('gzip', 5))

    plt.close(fig)
