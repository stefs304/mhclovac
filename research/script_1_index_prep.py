from pyaaisc import Aaindex
import joblib
from mhclovac.preprocessing import normalize_index_data, standardize_index_data


aaindex = Aaindex()
aaindex_data = {}

for accession_number, title in aaindex.get_all('aaindex1'):

    try:
        index_data = aaindex.get(accession_number, dbkey='aaindex1').index_data
        normalized_data = normalize_index_data(index_data)
        standardized_data = standardize_index_data(index_data)
        aaindex_data[accession_number] = {
            'title': title,
            'normalized_index_data': normalized_data,
            'standardized_index_data': standardized_data,
            'original_index_data': index_data
        }
        print(accession_number)
    except Exception as e:
        msg = f'{accession_number}: {e}'
        print(msg)


joblib.dump(aaindex_data, '../data/index_data.gz', compress=('gzip', 5))
