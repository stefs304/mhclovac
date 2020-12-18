import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

data = pd.read_csv('results/benchmarking_results.csv')

plt.figure(figsize=(10, 6), dpi=100)

sns.set_theme(style="whitegrid")
palette = {k: 'blue' for k in data.columns}
palette['mhclovac'] = 'red'
ax = sns.boxplot(data=data, color='#4c72b0')
ax.artists[0].set_facecolor('red')
ax.set_xticklabels([k for k in data.columns], rotation=30)
ax.set_yscale('log')
ax.set_ylim([10**-4, 1.0])
ax.set_ylabel('FRANK')
plt.savefig('figures/mhclovac-frank-benchmark.png', bbox_inches='tight')
