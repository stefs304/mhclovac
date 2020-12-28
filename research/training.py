from mhclovac.models import BindingModel
from mhclovac.preprocessing import get_features, transform_ic50_values
from mhclovac.config import Config
from sklearn.metrics import r2_score
import pandas as pd
import matplotlib.pyplot as plt


TRAINING_SET_FRACTION_LIST = [0.2, 0.4, 0.6, 0.8]
TRAINING_SET_FRACTION = 0.9
TRAINING_SET_SIZE_THRESHOLD = 50
RANDOM_SEED = 0
PLOT_PROPS = dict(boxstyle='round', facecolor='wheat', alpha=0.6)


data = pd.read_csv(f'../data/combined_data.zip')
for mhc_key in list(data['mhc_allele'].unique()):

    if mhc_key not in ['HLA-A*03:01']:
        continue

    mhc_data = data[data['mhc_allele'] == mhc_key]

    if len(mhc_data) < TRAINING_SET_SIZE_THRESHOLD:
        continue

    print(f'{mhc_key}')
    print(f'n_samples={len(mhc_data)}')

    fig, ax = plt.subplots(1, 3, constrained_layout=True, figsize=(9, 3))
    fig.suptitle(f'{mhc_key} training stats')

    # plot histogram of binding affinity measurements
    ax[0].hist(transform_ic50_values(mhc_data['affinity']), bins=50)
    # ax[0].set_xscale('log')
    ax[0].grid()
    ax[0].set_xlabel('log transformed affinity')
    ax[0].set_ylabel('n_samples')
    ax[0].text(0.05, 0.95, f'n_samples = {len(mhc_data)}', transform=ax[0].transAxes, verticalalignment='top', bbox=PLOT_PROPS)

    # calculate r2 score for different fractions of data used for training
    fraction_r2_scores = []
    for fraction in TRAINING_SET_FRACTION_LIST:

        train = mhc_data.sample(frac=fraction)
        test = mhc_data.drop(index=train.index)

        x_train = get_features(peptide_list=train['peptide'], index_id_list=Config.INDEX_ID_LIST)
        y_train = transform_ic50_values(train['affinity'])

        x_test = get_features(peptide_list=test['peptide'], index_id_list=Config.INDEX_ID_LIST)
        y_test = transform_ic50_values(test['affinity'])

        model = BindingModel(random_state=RANDOM_SEED)
        model.fit(x_train, y_train)

        predictions = model.predict(x_test)
        r2 = round(r2_score(y_test, predictions), 3)
        fraction_r2_scores.append(r2)

    ax[1].bar(TRAINING_SET_FRACTION_LIST, fraction_r2_scores, width=0.15)
    ax[1].grid()
    ax[1].set_ylim(0, 1)
    ax[1].set_xlabel('fraction of data used for training')
    ax[1].set_ylabel('r2_score')

    # train final model on whole dataset, preform testing on 0.1 of data
    train = mhc_data.sample(frac=1.0)
    test = mhc_data.sample(frac=0.1)

    x_train = get_features(peptide_list=train['peptide'], index_id_list=Config.INDEX_ID_LIST)
    y_train = transform_ic50_values(train['affinity'])

    x_test = get_features(peptide_list=test['peptide'], index_id_list=Config.INDEX_ID_LIST)
    y_test = transform_ic50_values(test['affinity'])

    model = BindingModel(random_state=RANDOM_SEED)
    model.fit(x_train, y_train)

    predictions = model.predict(x_test)
    r2 = round(r2_score(y_test, predictions), 3)
    print(f'r2_score={r2}')

    ax[2].scatter(y_test, predictions, s=4)
    ax[2].grid()
    ax[2].set_xlabel('log transformed affinity')
    ax[2].set_ylabel('predicted score')
    ax[2].text(0.05, 0.95, f'r2_score = {r2}', transform=ax[2].transAxes, verticalalignment='top', bbox=PLOT_PROPS)

    plt.savefig(f'./training_results/{mhc_key}.png', bbox_inches='tight')
    plt.close()



