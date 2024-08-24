
import pandas as pd
import numpy as np
import time
import math
import matplotlib.pyplot as plt
from mhclovac.binding import AlleleSpecificBindingPredictor



TRAINING_SET_SIZE_THRESHOLD = 50
RANDOM_SEED = 0

train_data = pd.read_csv('./data/mhc_full_cleaned.csv')

for mhc_allele in list(train_data['mhc_allele'].unique()):

    if not mhc_allele in ['HLA-B*27:05']:
        continue

    mhc_train_data = train_data[train_data['mhc_allele'] == mhc_allele]
    mhc_train_data.dropna(subset=['quant_meas'], inplace=True)
    target = mhc_train_data['quant_meas'].apply(lambda x: 1 - math.log(x, 50.000))

    if len(mhc_train_data) < TRAINING_SET_SIZE_THRESHOLD:
        continue

    print(f'training {mhc_allele}')

    model = AlleleSpecificBindingPredictor(mhc_allele=mhc_allele, median_length=9)
    model.train(mhc_train_data['peptide'], target)

    data = model.predict(mhc_train_data['peptide'])

    plt.scatter(target, data)
    plt.title(f'{np.corrcoef(target, data)}')
    plt.show()

