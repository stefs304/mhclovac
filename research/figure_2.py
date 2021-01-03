import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import joblib
from sklearn.metrics import auc


frank_data = pd.read_csv('results/benchmarking_frank.csv')
roc_data = joblib.load('results/benchmarking_roc.gzip')

sns.set_theme(style="whitegrid")
fig, ax = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)

sns.boxplot(y=frank_data['mhclovac'], whis=1.5, ax=ax[0], width=0.2)
ax[0].set_ylim(-0.05, 0.5)
ax[0].set_ylabel('FRANK')
ax[0].set_xlabel('MHCLovac 4.0')

auc_binding = round(auc(roc_data['fpr_array'], roc_data['tpr_array']), 3)
ax[1].plot(roc_data['fpr_array'], roc_data['tpr_array'], linewidth=3)
ax[1].set_xlabel('False positive rate')
ax[1].set_ylabel('True positive rate')
ax[1].legend(labels=(
    f'AUC = {auc_binding}',
))

plt.savefig('figures/mhclovac-benchmark.png', bbox_inches='tight')

