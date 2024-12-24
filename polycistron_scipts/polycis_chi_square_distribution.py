import pandas as pd
import numpy as np
import sys
from scipy.stats import chisquare
import matplotlib.pyplot as plt
bed_file = sys.argv[1]
df = pd.read_csv(bed_file, sep='\t', header=None, names=['chr', 'start', 'end'])
chr_sizes = {
    'CpBGF_Chr1': 919662,
    'CpBGF_Chr2': 991859,
    'CpBGF_Chr3': 1102334,
    'CpBGF_Chr4': 1107096,
    'CpBGF_Chr5': 1085428,
    'CpBGF_Chr6': 1308109,
    'CpBGF_Chr7': 1363375,
    'CpBGF_Chr8': 1379130
}


observed_distribution = df['chr'].value_counts().sort_index()


observed_distribution = observed_distribution[observed_distribution.index.isin(chr_sizes.keys())]

expected_distribution = pd.Series(chr_sizes).sort_index() / sum(chr_sizes.values()) * len(df)

expected_distribution = expected_distribution[expected_distribution.index.isin(chr_sizes.keys())]

chi2, p_value = chisquare(observed_distribution, expected_distribution)

print(f'Chi-square: {chi2}')
print(f'p-value: {p_value}')

fig, ax = plt.subplots()
bar_width = 0.35
index = np.arange(len(observed_distribution))

observed_bars = ax.bar(index, observed_distribution / len(df), bar_width, label='Observed')
expected_bars = ax.bar(index + bar_width, expected_distribution / len(df), bar_width, alpha=0.5, label='Expected')

ax.set_xlabel('Chromosome')
ax.set_ylabel('Proportion')
ax.set_xticks(index + bar_width / 2)
ax.set_xticklabels(observed_distribution.index, rotation=45, ha='right')  # Rotate x-axis labels at a 45-degree angle
ax.legend()

def autolabel(bars):
    for bar in bars:
        height = bar.get_height()
        ax.annotate('{:.2f}'.format(height),
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

autolabel(observed_bars)
autolabel(expected_bars)

ax.annotate('p-value: {:.4f}\nChi-square: {:.4f}'.format(p_value, chi2),
            xy=(0.5, 0.95), xycoords='axes fraction',
            xytext=(0, -30), textcoords='offset points',
            ha='center', va='top')

fig.savefig('chi_square_plots.png', dpi=600)
plt.close(fig)
