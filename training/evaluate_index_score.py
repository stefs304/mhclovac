import json
import os
import pickle
from mhclovac.utils import index_correlation
from pprint import pprint


training_configs = {}
for fname in os.listdir('./exploration_results'):
    with open(os.path.join('exploration_results', fname)) as f:
        mhc_allele = fname.split('.json')[0]
        config = json.load(f)
        if not (config['binding'] or config['epitope']):
            continue
        if mhc_allele == config['mhc_key']:
            training_configs[mhc_allele] = config
        else:
            print('file name doesn\'t match mhc_key')


with open('../data/aaindex_data.pickle', 'rb') as f:
    aaindex_data = pickle.load(f)


sum_score = {}
count = 0
for key in training_configs:
    if training_configs[key]['binding']:
        count += 1
        for entry in training_configs[key]['binding']:
            if entry[0] in sum_score:
                sum_score[entry[0]] += entry[2]
            else:
                sum_score[entry[0]] = entry[2]


sum_score = [(key, val / count) for key, val in sum_score.items()]
sum_score = sorted(sum_score, key=lambda x: x[1], reverse=True)

correlation_threshold = 0.3
index_list = []
count = 0
for new_entry in sum_score:
    # if count >= 10:
        # break
    if not index_list:
        index_list.append(new_entry)
        count += 1
        continue
    flag = True
    for entry in index_list:
        correlation = index_correlation(
            aaindex_data[new_entry[0]]['index_data'],
            aaindex_data[entry[0]]['index_data']
        )
        if correlation > correlation_threshold or correlation < -correlation_threshold:
            flag = False
    if flag:
        if new_entry[1] > 0:
            index_list.append(new_entry)
            count += 1


details = []
for entry in index_list:
    if entry[1] > 0:
        details.append((entry[0], aaindex_data[entry[0]]['title'], entry[1]))
        print(entry[0], aaindex_data[entry[0]]['title'], entry[1])

print('\n')
print([x[0] for x in index_list])

id_list = [x[0] for x in index_list]

print(id_list)
