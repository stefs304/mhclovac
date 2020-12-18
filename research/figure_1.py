import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
from mhclovac.utils import load_index_data, pdf
from mhclovac.sequence import model_distribution
import seaborn as sns

sns.set_theme(style="whitegrid")
# Setup plot formatting, fonts etc.
rcParams['font.family'] = 'sans-serif'
# rcParams['font.sans-serif'] = ['Tahoma']
rcParams['font.size'] = 16
props = dict(boxstyle='round', facecolor='wheat', alpha=0.6)


hydrophobicity_index_data = load_index_data()['HOPT810101']['standardized_index_data']


def plot_distribution(sequence, index, label, legend=None, linewidth=3, ylabel=None, title=None):
    multiplier = 24
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
        plt.plot(current_vector, '--', linewidth=linewidth)
    k, = plt.plot(dist_vector, 'k', linewidth=linewidth)
    plt.xticks(x_ticks, x_labels)
    if ylabel:
        plt.ylabel(ylabel, fontsize='large')
    y = max(dist_vector) - max(dist_vector) / 20.0
    plt.text(0, y, label, fontsize='large')
    if legend:
        plt.legend([k],[legend])
    if title:
        plt.title(title)
    return


def plot_discrete_values(distribution, title=None, label=None):
    n_discrete_points = 10
    step = int(len(distribution) / n_discrete_points)
    x_axis = []
    y_axis = []
    for i in range(0, len(distribution), step):
        x_axis.append(i)
        y_axis.append(distribution[i])

    plt.bar(range(n_discrete_points), y_axis)
    plt.xticks(range(n_discrete_points), range(n_discrete_points))
    y = max(y_axis) - max(y_axis) / 20.0
    if label:
        plt.text(0, y, label, fontsize='large')
    if title:
        plt.title(title)
    return


example_sequence_1 = 'LLDVTAAV'
example_sequence_2 = 'FLFDGSPTYVL'


plt.figure(figsize=(12, 12), dpi=100)

plt.subplot(3, 1, 1)
plot_distribution(
    example_sequence_1,
    hydrophobicity_index_data,
    label='a)',
    title='Modeled hydrophilicity profile'
)

dist_1 = model_distribution(example_sequence_1, hydrophobicity_index_data)

plt.subplot(3, 2, 3)
plot_distribution(
    example_sequence_1,
    hydrophobicity_index_data,
    label='b)',
    # legend=example_sequence_1,
    linewidth=2,
    title='8-mer profile',
    ylabel='Hydrophilicity index'
)

plt.subplot(3, 2, 5)
plot_discrete_values(
    distribution=dist_1,
    title='8-mer discrete profile',
    label='d)'
)


dist_2 = model_distribution(example_sequence_2, hydrophobicity_index_data)
plt.subplot(3, 2, 4)
plot_distribution(
    example_sequence_2,
    hydrophobicity_index_data,
    label='c)',
    # legend=example_sequence_2,
    linewidth=2,
    title='11-mer profile'
)

plt.subplot(3, 2, 6)
plot_discrete_values(
    distribution=dist_2,
    label='e)',
    title='11-mer discrete profile'
)


plt.savefig('figures/mhclovac-modeling-figure.png', bbox_inches='tight')
