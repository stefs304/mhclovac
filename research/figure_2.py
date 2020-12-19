import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
import joblib
from sklearn.metrics import auc


frank_data = pd.read_csv('results/benchmarking_frank.csv')
roc_data = joblib.load('results/benchmarking_roc.gzip')

sns.set_theme(style="whitegrid")

fig = plt.figure(constrained_layout=True)
gs = GridSpec(1, 3, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1:])

sns.boxplot(y=frank_data['mhclovac'], whis=1.5, ax=ax1)
# sns.swarmplot(data=frank_data, color=".25", ax=ax1)

# ax1.set_xticklabels(['MHCLovac 4.0'])
# ax1.set_yscale('log')
ax1.set_ylim(-0.05, 0.5)
ax1.set_ylabel('FRANK')
ax1.set_xlabel('MHCLovac 4.0')

auc_binding = round(auc(roc_data['fpr_array'], roc_data['tpr_array']), 3)
ax2.plot(roc_data['fpr_array'], roc_data['tpr_array'], linewidth=3)
ax2.set_xlabel('False positive rate')
ax2.set_ylabel('True positive rate')
ax2.legend(labels=(
    f'AUC = {auc_binding}',
))


plt.savefig('figures/mhclovac-benchmark.png', bbox_inches='tight')

