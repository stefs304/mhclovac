import pandas as pd
import numpy as np
import pickle
import json
from mhclovac.preprocessing import sequence_to_features
from mhclovac.models import BindingPredictor, EpitopePredictor
from sklearn.metrics import r2_score, f1_score
import logging
import multiprocessing as mp


logging.basicConfig(filename='exp.log', level=logging.DEBUG, format='%(asctime)s: %(message)s')
logging.info(f'Starting exploration')


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


train_set_percentage = 0.9

train_set_size_threshold = 50

percentile_threshold = 95

n_epochs = 3

n_processes = 2

skip_list = [
    'HLA class I',
    'B19 class I',
    'H2-Kb Y22F, M23I, E24S, D30N mutant',
    'RT1-A',
]

with open('../data/aaindex_data.pickle', 'rb') as f:
    aaindex_data = pickle.load(f)


data = pd.read_csv(f'../data/combined_data.tsv', sep='\t', low_memory=False)
data['meas'] = data['meas'].apply(lambda x: 20000 if x > 20000 else x)  # cap at 20k

# -----------------------------------------------------------------------------------------


def binding_worker(index_id, data, aaindex_data, n_epochs, result_queue):
    index_prediction_scores = []
    index_data = [aaindex_data[index_id]['index_data']]
    for _ in range(n_epochs):

        train = data.sample(frac=train_set_percentage)
        test = data.drop(index=train.index)

        X_train = train['sequence'].apply(lambda x: sequence_to_features(x, index_data))
        X_train = pd.DataFrame(X_train.tolist(), index=train.index)
        y_train = train['meas'].apply(meas_func)

        try:
            model = BindingPredictor(random_state=None, standardize_data=True, n_jobs=1)
            model.fit(X_train, y_train)

            X_test = test['sequence'].apply(lambda x: sequence_to_features(x, index_data))
            X_test = pd.DataFrame(X_test.tolist(), index=test.index)
            y_test = test['meas'].apply(meas_func)

            predictions = model.predict(X_test)
            r2 = round(r2_score(y_test, predictions), 3)

            index_prediction_scores.append(r2)

        except Exception as e:
            logging.warning(e)
            pass

    mean_prediction_score = np.mean(index_prediction_scores)
    result_data = (index_id, aaindex_data[index_id]['title'], mean_prediction_score)
    result_queue.put(result_data)
    return None


def epitope_worker(index_id, data, aaindex_data, n_epochs, result_queue):
    index_prediction_scores = []
    index_data = [aaindex_data[index_id]['index_data']]
    for _ in range(n_epochs):

        train = data.sample(frac=train_set_percentage)
        test = data.drop(index=train.index)

        X_train = train['sequence'].apply(lambda x: sequence_to_features(x, index_data))
        X_train = pd.DataFrame(X_train.tolist(), index=train.index)
        y_train = train['label'].apply(label_func)

        try:
            model = EpitopePredictor(random_state=None, standardize_data=True, n_jobs=1)
            model.fit(X_train, y_train)

            X_test = test['sequence'].apply(lambda x: sequence_to_features(x, index_data))
            X_test = pd.DataFrame(X_test.tolist(), index=test.index)
            y_test = test['label'].apply(label_func)

            preds_binary = model.predict(X_test)
            f1 = round(f1_score(y_test, preds_binary), 3)

            index_prediction_scores.append(f1)

        except Exception as e:
            logging.warning(e)
            pass

    mean_prediction_score = np.mean(index_prediction_scores)
    result_data = (index_id, aaindex_data[index_id]['title'], mean_prediction_score)
    result_queue.put(result_data)
    return None


for mhc_key in list(data['mhc'].unique()):

    # if mhc_key != 'Patr-B*13:01':
    #    continue

    if mhc_key in skip_list:
        logging.info(f'skipping {mhc_key}')
        continue

    logging.info(mhc_key)

    exploration_results = {
        'mhc_key': mhc_key,
        'binding': None,
        'epitope': None
    }

    tmp_data = data[data['mhc'] == mhc_key]
    logging.info(f'Data shape: {tmp_data.shape}')

    tmp_bdata = tmp_data[['sequence', 'meas']].dropna().copy()

    if len(tmp_bdata) >= train_set_size_threshold:

        logging.info(f'exploring binding prediction')

        manager = mp.Manager()
        result_queue = manager.Queue()
        pool = mp.Pool(processes=n_processes)

        for i, index_id in enumerate(aaindex_data):
            r = pool.apply_async(binding_worker, args=(index_id, tmp_bdata, aaindex_data, n_epochs, result_queue,), error_callback=print)

        pool.close()
        pool.join()

        result_queue.put('DONE')

        result_list = []
        while True:
            result_data = result_queue.get()
            if result_data == 'DONE':
                break
            result_list.append(result_data)

        exploration_results['binding'] = result_list

        # sad treba naci 95 percentil i pronaci one indekse koji prolaze

        # tmp_perc_threshold = np.percentile([x for x in prediction_results.values()], percentile_threshold)
        # pass_list = []
        #
        # for index_id, mean_score in prediction_results.items():
        #     if mean_score >= tmp_perc_threshold:
        #         pass_list.append([index_id, aaindex_data[index_id]['title'], mean_score])
        #
        # exploration_results['binding'] = pass_list

    # sad isto za epitope
    tmp_edata = tmp_data[['sequence', 'label']].dropna().copy()

    if len(tmp_edata) >= train_set_size_threshold:

        logging.info(f'exploring epitope prediction')

        manager = mp.Manager()
        result_queue = manager.Queue()
        pool = mp.Pool(processes=n_processes)

        for i, index_id in enumerate(aaindex_data):
            pool.apply_async(epitope_worker, args=(index_id, tmp_edata, aaindex_data, n_epochs, result_queue,), error_callback=print)

        pool.close()
        pool.join()

        result_queue.put('DONE')

        result_list = []
        while True:
            result_data = result_queue.get()
            if result_data == 'DONE':
                break
            result_list.append(result_data)

        exploration_results['epitope'] = result_list

        # sad treba naci 95 percentil i pronaci one indekse koji prolaze

        # tmp_perc_threshold = np.percentile([x for x in prediction_results.values()], percentile_threshold)
        # pass_list = []
        #
        # for index_id, mean_score in prediction_results.items():
        #     if mean_score >= tmp_perc_threshold:
        #         pass_list.append([index_id, aaindex_data[index_id]['title'], mean_score])
        #
        # exploration_results['epitope'] = pass_list

    with open(f'exploration_results/{mhc_key}.json', 'w') as f:
        f.write(json.dumps(exploration_results, indent=2))
        logging.info(f'file "{mhc_key}.json" saved')

