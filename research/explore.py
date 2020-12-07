from mhclovac.models import BindingModel
from mhclovac.preprocessing import get_features, transform_ic50_values
from mhclovac.utils import load_index_data, get_index_correlation
from sklearn.metrics import r2_score
import numpy as np
import pandas as pd
import multiprocessing as mp
import logging
import joblib


TRAINING_SET_FRACTION = 0.8
TRAINING_SET_SIZE_THRESHOLD = 50
DATA_FILE = '../data/combined_data.tsv'
SKIP_LIST = [
    'HLA class I',
    'B19 class I',
    'H2-Kb Y22F, M23I, E24S, D30N mutant',
    'RT1-A',
]
N_PROC = 6
N_ITERATIONS = 5
CORRELATION_THRESHOLD = 0.3


logging.basicConfig(filename='explore.log', level=logging.DEBUG, format='%(asctime)s: %(message)s')

data = pd.read_csv(DATA_FILE, sep='\t', low_memory=False)
index_data = load_index_data()


def worker(data, index_id, n_iterations, result_queue):
    index_prediction_scores = []
    for _ in range(n_iterations):
        train = data.sample(frac=TRAINING_SET_FRACTION)
        test = data.drop(index=train.index)

        x_train = get_features(peptide_list=train['sequence'], index_id_list=[index_id])
        y_train = transform_ic50_values(train['meas'])

        try:
            model = BindingModel(n_jobs=1)
            model.fit(x_train, y_train)

            x_test = get_features(peptide_list=test['sequence'], index_id_list=[index_id])
            y_test = transform_ic50_values(test['meas'])

            predictions = model.predict(x_test)
            r2 = round(r2_score(y_test, predictions), 3)

            index_prediction_scores.append(r2)

        except Exception as e:
            logging.error(e)
            pass

    mean_prediction_score = np.mean(index_prediction_scores)
    result_data = (index_id, mean_prediction_score)
    result_queue.put(result_data)
    return None


exploration_output = {}

for mhc_key in list(data['mhc'].unique()):

    if mhc_key not in ['Patr-B*13:01', 'H2-Kb']:
        continue

    if mhc_key in SKIP_LIST:
        logging.info(f'skipping {mhc_key}')
        continue

    logging.info(mhc_key)

    mhc_results = {
        'mhc_key': mhc_key,
        'results': None
    }

    tmp_data = data[data['mhc'] == mhc_key]
    tmp_bdata = tmp_data[['sequence', 'meas']].dropna().copy()

    if len(tmp_bdata) < TRAINING_SET_SIZE_THRESHOLD:
        continue

    logging.info(f'Data shape: {tmp_data.shape}')
    manager = mp.Manager()
    result_queue = manager.Queue()
    pool = mp.Pool(processes=N_PROC)

    for i, index_id in enumerate(index_data):
        if i > 5: break
        r = pool.apply_async(worker, args=(tmp_bdata, index_id, N_ITERATIONS, result_queue,), error_callback=print)

    pool.close()
    pool.join()

    result_queue.put('DONE')

    result_list = []
    while True:
        result_data = result_queue.get()
        if result_data == 'DONE':
            break
        result_list.append(result_data)

    mhc_results['results'] = result_list
    exploration_output[mhc_key] = mhc_results


joblib.dump(exploration_output, 'exploration_output.gz', compress=('gzip', 5))
logging.info(f'Data saved.')

# --------------

index_scores = {}
for index_id in index_data:

    for mhc_key in exploration_output:
        for prediction_result in exploration_output[mhc_key]['results']:
            if prediction_result[0] in index_scores:
                index_scores[prediction_result[0]].append(prediction_result[1])
            else:
                index_scores[prediction_result[0]] = [prediction_result[1]]


average_scores = []
for index_id in index_scores:
    mean_score = np.mean(index_scores[index_id])
    average_scores.append((index_id, mean_score))

sorted_average_scores = sorted(average_scores, key=lambda x: x[1], reverse=True)

pass_index_list = []
for entry in sorted_average_scores:
    if not pass_index_list:
        pass_index_list.append(entry)
        continue
    flag = True
    for passed_entry in pass_index_list:
        correlation = get_index_correlation(
            index_data[entry[0]]['index_data'],
            index_data[passed_entry[0]]['index_data']
        )
        if np.abs(correlation) > CORRELATION_THRESHOLD:
            flag = False
    if flag and (entry[1] > 0):
        pass_index_list.append(entry)

logging.info(f'Pass index IDs')

for index_id, average_prediction_score in pass_index_list:
    title = index_data[index_id]['title']
    logging.info(f'{index_id} - {title} - {average_prediction_score}')

joblib.dump(pass_index_list, 'pass_index_list.gz', compress=('gzip', 5))

