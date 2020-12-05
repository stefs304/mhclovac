from pyaaisc import Aaindex
import pickle
from mhclovac.preprocessing import standardize_index


aaindex = Aaindex()

aaindex_data = {}

for accession_number, title in aaindex.get_all('aaindex1'):

    try:
        index_data = aaindex.get(accession_number, dbkey='aaindex1').index_data
        standardized_data = standardize_index(index_data)
        aaindex_data[accession_number] = {
            'title': title,
            'index_data': standardized_data
        }
    except Exception as e:
        msg = f'{accession_number}: {e}'
        print(msg)

with open('../data/aaindex_data.pickle', 'wb') as f:
    pickle.dump(aaindex_data, f)
