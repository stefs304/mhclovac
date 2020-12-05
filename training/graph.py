import pickle
import pandas as pd
import numpy as np
from proteinko import encode_sequence
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from mhclovac.models import BindingPredictor


with open('../data/aaindex_data.pickle', 'rb') as f:
    aaindex_data = pickle.load(f)


def meas_func(x):
    return (10 - np.log(x)) / 10


def normalize_index(index_data):
    values = [val for val in index_data.values()]
    for key, val in index_data.items():
        new_val = (val - min(values)) / (max(values) - min(values))
        new_val = 2 * new_val - 1
        index_data[key] = new_val
    return index_data


for key in aaindex_data:
    normalize_index(aaindex_data[key]['index_data'])


def get_bin(value):
    ranges = [
        (-1.0, -0.66),
        (-0.66, -0.33),
        (-0.33, 0),
        (0, 0.33),
        (0.33, 0.66),
        (0.66, 1.1)     # 1.1 for logic
    ]
    for range_ in ranges:
        if value >= range_[0] and value < range_[1]:
            return str(range_)
    return None


index_to_use = ['PRAM900101', 'FASG760101', 'ZIMJ680104', 'CHOP780201', 'PRAM820103', 'RACS820112', 'ROBB760107']


def sequence_to_graph(sequence, score):
    node_list = []
    graph = {}
    for index in index_to_use:
        index_data = aaindex_data[index]['index_data']
        encoded_seq = encode_sequence(sequence=sequence, encoding_scheme=index_data)
        for i, value in enumerate(encoded_seq):
            value_bin = get_bin(value)
            node_name = f'{index}:{i}:{value_bin}'
            node_list.append(node_name)
    for i in range(0, len(node_list)-1):
        for j in range(i+1, len(node_list)):
            edge = (node_list[i], node_list[j])
            graph[edge] = score
    return graph


def sequence_to_graph_array(sequence_graph, reference_graph):
    array = []
    for key in reference_graph:
        if key in sequence_graph:
            array.append(1)
        else:
            array.append(0)
    return array


def predict_score(sequence_graph, reference_graph):
    score = 0
    for edge in sequence_graph:
        if edge in reference_graph:
            score += reference_graph[edge]
    return score


data = pd.read_csv(f'../data/combined_data.tsv', sep='\t')
data['meas'] = data['meas'].apply(lambda x: 20000 if x > 20000 else x)  # cap at 20k
data['meas'] = data['meas'].apply(lambda x: 1.0 if x < 1.0 else x)


for mhc_key in list(data['mhc'].unique()):
    if mhc_key != 'HLA-A*24:02':
        continue

    tmp_data = data[data['mhc'] == mhc_key]
    tmp_data = tmp_data[['sequence', 'meas']].dropna().copy()

    tmp_data['score'] = tmp_data['meas'].apply(meas_func)
    tmp_data['weights'] = tmp_data['meas'].apply(lambda x: 1/x)
    tmp_data['sequence_length'] = tmp_data['sequence'].apply(len)

    train = tmp_data.sample(frac=0.9)
    test = tmp_data.drop(index=train.index)

    print('building graph')
    graph = {}
    for i, row in train.iterrows():
        row_graph = sequence_to_graph(row['sequence'], row['weights'])
        for edge in row_graph:
            if edge in graph:
                graph[edge] += row_graph[edge]
            else:
                graph[edge] = row_graph[edge]
    for key in graph:
        graph[key] = graph[key] / float(len(tmp_data))
    print(f'graph size: {len(graph)}')

    print('sorting graph')
    sorted_graph = {k: v for k, v in sorted(graph.items(), key=lambda item: item[1], reverse=True)}
    top_graph = {}
    for i, key in enumerate(sorted_graph):
        if i > 100:
            break
        top_graph[key] = sorted_graph[key]

    # predictions = []
    # for i, row in tmp_data.iterrows():
    #     seq_graph = sequence_to_graph(row['sequence'], 0)
    #     predictions.append(predict_score(seq_graph, graph))
    # # print(predictions)
    # r2 = r2_score(tmp_data['score'], predictions)
    # print(f'r2: {r2}')
    # plt.scatter(tmp_data['score'], predictions)
    # plt.savefig('test.png')

    features = []
    for i, row in train.iterrows():
        seq_graph = sequence_to_graph(row['sequence'], 0)
        feat_array = sequence_to_graph_array(seq_graph, top_graph)
        features.append(feat_array)

    features = pd.DataFrame(features)
    model = BindingPredictor(standardize_data=False)
    model.fit(features, train['score'])

    test_features = []
    for i, row in test.iterrows():
        seq_graph = sequence_to_graph(row['sequence'], 0)
        feat_array = sequence_to_graph_array(seq_graph, top_graph)
        test_features.append(feat_array)

    predictions = model.predict(test_features)

    r2 = r2_score(test['score'], predictions)
    print(f'r2: {r2}')

    plt.figure(figsize=(5, 5))
    plt.scatter(test['score'], predictions, s=4)
    plt.grid()
    plt.savefig('test.png')
