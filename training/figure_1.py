import matplotlib.pyplot as plt
import numpy as np
from mhclovac.utils import load_index_data, pdf
from mhclovac.sequence import model_distribution
import seaborn as sns


def plot_distribution(ax, sequence, index, ylabel=None, title=None):
    multiplier = 20
    overlap_distance = 1
    sigma = 0.8
    dist_vector = np.zeros(multiplier * len(sequence) + 2 * overlap_distance * multiplier)
    x_ticks = []
    x_labels = []
    for i, aa in enumerate(sequence):
        current_vector = np.zeros(multiplier * len(sequence) + 2 * overlap_distance * multiplier)
        value = index[aa]
        x = np.linspace(-2.3263, 2.3263, (2 * overlap_distance + 1) * multiplier)
        aa_dist = pdf(x, sigma) * value
        current_vector[int(i * multiplier):int((i + (2 * overlap_distance + 1)) * multiplier)] = aa_dist
        dist_vector[int(i * multiplier):int((i + (2 * overlap_distance + 1)) * multiplier)] += aa_dist
        x_ticks.append(np.argmax(np.abs(current_vector)))
        x_labels.append(aa)
        ax.plot(current_vector, '--', linewidth=2)
    ax.plot(dist_vector, 'k', linewidth=3)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    return


def plot_discrete_profile(ax, sequence, index, title=None, ylabel=None):
    discrete_profile = model_distribution(sequence, index, n_discrete_points=10)
    ax.bar(range(len(discrete_profile)), discrete_profile)
    ax.set_xticks(range(len(discrete_profile)))
    ax.set_xticklabels(range(len(discrete_profile)))
    if title:
        ax.set_title(title)
    if ylabel:
        ax.set_ylabel(ylabel)
    return


sns.set_theme(style="whitegrid")
hydrophobicity_index_data = load_index_data()['ROSM880102']['standardized_index_data']
example_sequence_1 = 'LLDVTAAV'
example_sequence_2 = 'FLFDGSPTYVL'


fig, ax = plt.subplots(2, 2, figsize=(10, 5), constrained_layout=True)

plot_distribution(
    ax=ax[0][0],
    sequence='LLDVTAAV',
    index=hydrophobicity_index_data,
    title='8-mer profile',
    ylabel='Hydropathy index'
)

plot_distribution(
    ax=ax[0][1],
    sequence='FLFDGSPTYVL',
    index=hydrophobicity_index_data,
    title='11-mer profile'
)

plot_discrete_profile(
    ax=ax[1][0],
    sequence='LLDVTAAV',
    index=hydrophobicity_index_data,
    title='8-mer discrete profile',
    ylabel='Hydropathy index'
)

plot_discrete_profile(
    ax=ax[1][1],
    sequence='FLFDGSPTYVL',
    index=hydrophobicity_index_data,
    title='11-mer discrete profile'
)


plt.savefig('figures/mhclovac-modeling-figure.png', bbox_inches='tight')
