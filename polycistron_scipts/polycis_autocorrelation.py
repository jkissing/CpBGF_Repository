import sys
import numpy as np
import matplotlib.pyplot as plt

bed_file = sys.argv[1]

features = []
with open(bed_file, 'r') as f:
    for line in f:
        chrom, start, end = line.strip().split('\t')
        features.append((chrom, int(start), int(end)))

chrom_sizes = {
    'CpBGF_Chr1': 919662,
    'CpBGF_Chr2': 991859,
    'CpBGF_Chr3': 1102334,
    'CpBGF_Chr4': 1107096,
    'CpBGF_Chr5': 1085428,
    'CpBGF_Chr6': 1308109,
    'CpBGF_Chr7': 1363375,
    'CpBGF_Chr8': 1379130
}

max_distance_threshold_fraction = 0.1

num_thresholds = 10

num_bins = 20

max_distance_threshold = min(chrom_sizes.values()) * max_distance_threshold_fraction

distances_by_chrom = {}
for chrom in chrom_sizes.keys():
    chrom_features = [feature for feature in features if feature[0] == chrom]
    distances = []
    for feature1 in chrom_features:
        start1, _ = feature1[1:]
        for feature2 in chrom_features:
            start2, _ = feature2[1:]
            distance = abs(start1 - start2)
            if distance <= max_distance_threshold:
                distances.append(distance)
    distances_by_chrom[chrom] = distances

num_simulations = 1000
correlations_observed_by_chrom = {}
correlations_simulated_by_chrom = {}
p_values_by_chrom = {}

for chrom in chrom_sizes.keys():
    distances = distances_by_chrom[chrom]

    best_threshold = None
    best_bin_number = None
    best_p_value = 1.0

    for i in range(1, num_thresholds + 1):
        distance_threshold = max_distance_threshold * (i / num_thresholds)

        max_distance = max(distances)
        bin_size = max_distance / num_bins
        bins = np.arange(0, max_distance + bin_size, bin_size)
        correlations_observed = np.zeros(len(bins) - 1)
        for j, bin_start in enumerate(bins[:-1]):
            bin_end = bin_start + bin_size
            counts = sum(bin_start <= distance <= bin_end for distance in distances)
            correlations_observed[j] = counts

        correlations_simulated = np.zeros((num_simulations, len(correlations_observed)))
        for k in range(num_simulations):
            np.random.shuffle(distances)
            correlations = np.zeros(len(bins) - 1)
            for j, bin_start in enumerate(bins[:-1]):
                bin_end = bin_start + bin_size
                counts = sum(bin_start <= distance <= bin_end for distance in distances)
                correlations[j] = counts
            correlations_simulated[k] = correlations

        p_values = np.mean(correlations_simulated >= correlations_observed, axis=0)
        min_p_value = np.min(p_values)
        if min_p_value < best_p_value:
            best_threshold = distance_threshold
            best_bin_number = num_bins
            best_p_value = min_p_value

    p_values_by_chrom[chrom] = p_values
    correlations_observed_by_chrom[chrom] = correlations_observed
    correlations_simulated_by_chrom[chrom] = correlations_simulated

fig, axs = plt.subplots(nrows=len(chrom_sizes), ncols=2, figsize=(12, 6 * len(chrom_sizes)), sharex=True)

for i, chrom in enumerate(chrom_sizes.keys()):
    distances = distances_by_chrom[chrom]
    correlations_observed = correlations_observed_by_chrom[chrom]
    correlations_simulated = correlations_simulated_by_chrom[chrom]
    p_values = p_values_by_chrom[chrom]

    bins = np.arange(0, len(correlations_observed))
    ax1 = axs[i, 0]
    ax2 = axs[i, 1]

    ax1.plot(bins, correlations_observed, marker='o', label='Observed')
    ax1.set_ylabel('Correlation Coefficient')
    ax1.set_title(f'Spatial Autocorrelation - Chromosome {chrom}')

    for x, y, p in zip(bins, correlations_observed, p_values):
        if p <= 0.01:
            ax1.text(x, y, f'{p:.2e}', fontsize=10, va='bottom', ha='center')
        elif p <= 0.05:
            ax1.text(x, y, f'{p:.2e}', fontsize=6, va='bottom', ha='center', color='gray')

    ax2.hist(correlations_simulated.flatten(), bins=30, alpha=0.5, label='Expected')
    ax2.axvline(correlations_observed.mean(), color='r', linestyle='dashed', linewidth=2, label='Observed')
    ax2.set_xlabel('Correlation Coefficient')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Expected vs Observed Correlations')
    ax2.legend()

    fig.tight_layout(pad=2.0)

fig.savefig('spatial_autocorrelation_plots.png', dpi=600)
plt.close(fig)
