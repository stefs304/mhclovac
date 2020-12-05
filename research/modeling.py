import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
from mhclovac.utils import load_index_data, pdf
from mhclovac.sequence import model_distribution

# Setup plot formatting, fonts etc.
rcParams['font.family'] = 'sans-serif'
# rcParams['font.sans-serif'] = ['Tahoma']
rcParams['font.size'] = 14
props = dict(boxstyle='round', facecolor='wheat', alpha=0.6)


hydrophobicity_index_data = load_index_data()['PRAM900101']['index_data']


def plot_distribution(sequence, index, label, legend=None, linewidth=3, ylabel=None):
    multiplier = 40
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
    plt.grid()
    if ylabel:
        plt.ylabel(ylabel)
    y = max(dist_vector) - max(dist_vector) / 20.0
    plt.text(0, y, label)
    if legend:
        plt.legend([k],[legend])


example_sequence_1 = 'LLDVTAAV'
example_sequence_2 = 'FLFDGSPTYVL'


plt.figure(figsize=(10, 10))

plt.subplot(3, 1, 1)
plot_distribution(example_sequence_1, hydrophobicity_index_data, label='a)', ylabel='Hydrophobicity')

n_discrete_points = 10
dist_1 = model_distribution(example_sequence_1, hydrophobicity_index_data)

plt.subplot(3, 2, 3)
plot_distribution(example_sequence_1, hydrophobicity_index_data, label='b)', legend=example_sequence_1, linewidth=2)

step = int(len(dist_1) / n_discrete_points)
x_axis = []
y_axis = []
for i in range(0, len(dist_1), step):
    x_axis.append(i)
    y_axis.append(dist_1[i])

plt.subplot(3, 2, 5)
plt.bar(range(10), y_axis)
plt.xticks(range(10), range(10))
plt.grid()
plt.text(0, 0.8, 'd)')


dist_2 = model_distribution(example_sequence_2, hydrophobicity_index_data)
plt.subplot(3, 2, 4)
plot_distribution(example_sequence_2, hydrophobicity_index_data, label='c)', legend=example_sequence_2, linewidth=2)


step = int(len(dist_2) / n_discrete_points)
x_axis = []
y_axis = []
for i in range(0, len(dist_2), step):
    x_axis.append(i)
    y_axis.append(dist_2[i])

plt.subplot(3, 2, 6)
plt.bar(range(10), y_axis)
plt.xticks(range(10), range(10))
plt.grid()
plt.text(0, 0.8, 'e)')


plt.savefig('figures/mhclovac-modeling-figure.png', bbox_inches='tight')
