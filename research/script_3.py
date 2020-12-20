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
DATA_FILE = '../data/combined_new.zip'
N_PROC = 2
N_ITERATIONS = 5
CORRELATION_THRESHOLD = 0.3

#create a logger
logger = logging.getLogger('exp')
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler('exp.log')
logger.addHandler(handler)


data = pd.read_csv(DATA_FILE)
index_data = load_index_data()


def worker(data, index_id, n_iterations, result_queue):
    logger = logging.getLogger('exp')
    index_prediction_scores = []
    mhc_key = list(data['mhc'].unique())[0]

    for _ in range(n_iterations):
        train = data.sample(frac=TRAINING_SET_FRACTION)
        test = data.drop(index=train.index)

        x_train = get_features(peptide_list=train['peptide'], index_id_list=[index_id])
        # y_train = transform_ic50_values(train['ic50'])
        y_train = train['ic50']
        try:
            model = BindingModel(n_jobs=1)
            model.fit(x_train, y_train)

            x_test = get_features(peptide_list=test['peptide'], index_id_list=[index_id])
            # y_test = transform_ic50_values(test['ic50'])
            y_test = test['ic50']

            predictions = model.predict(x_test)
            r2 = round(r2_score(y_test, predictions), 3)

            index_prediction_scores.append(r2)

        except Exception as e:
            logging.error(e)
            pass

    mean_prediction_score = np.mean(index_prediction_scores)
    result_data = (index_id, mean_prediction_score)
    logger.debug(f'{mhc_key}:{index_id}:{mean_prediction_score}')
    result_queue.put(result_data)
    return None


exploration_output = {}

for mhc_key in list(data['mhc'].unique()):

    logger.info(mhc_key)

    mhc_results = {
        'mhc_key': mhc_key,
        'results': None
    }

    mhc_data = data[data['mhc'] == mhc_key]

    if len(mhc_data) < TRAINING_SET_SIZE_THRESHOLD:
        continue

    logger.info(f'Data shape: {mhc_data.shape}')
    manager = mp.Manager()
    result_queue = manager.Queue()
    pool = mp.Pool(processes=N_PROC)

    for i, index_id in enumerate(index_data):
        r = pool.apply_async(worker, args=(mhc_data, index_id, N_ITERATIONS, result_queue,), error_callback=print)

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


joblib.dump(exploration_output, 'results/exploration_output.gz', compress=('gzip', 5))
logger.info(f'Data saved.')

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
            index_data[entry[0]]['standardized_index_data'],
            index_data[passed_entry[0]]['standardized_index_data']
        )
        if np.abs(correlation) > CORRELATION_THRESHOLD:
            flag = False
    if flag and (entry[1] > 0):
        pass_index_list.append(entry)

logger.info(f'Pass index IDs')

for index_id, average_prediction_score in pass_index_list:
    title = index_data[index_id]['title']
    logger.info(f'{index_id} - {title} - {average_prediction_score}')

joblib.dump(pass_index_list, 'results/pass_index_list.gz', compress=('gzip', 5))

