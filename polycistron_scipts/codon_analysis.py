import pandas as pd

def read_fasta(file_path):
    """Manually parse a FASTA file and return a list of sequences."""
    sequences = []
    current_sequence = ''
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_sequence:
                    sequences.append(current_sequence)
                    current_sequence = ''
            else:
                current_sequence += line
        if current_sequence:
            sequences.append(current_sequence)
    return sequences

def count_codons(sequences):
    """Count the codons in a list of sequences."""
    codon_counts = {}
    total_codons = 0
    for sequence in sequences:
        for i in range(0, len(sequence), 3):
            codon = sequence[i:i + 3].upper()
            if len(codon) == 3:
                codon_counts[codon] = codon_counts.get(codon, 0) + 1
                total_codons += 1
    return codon_counts, total_codons

def calculate_percentage(codon_counts, total_codons):
    """Calculate the percentage of each codon."""
    return {codon: (count / total_codons) * 100 for codon, count in codon_counts.items()}

polycistron_path = 'polycistron_cds.fasta'
emc_path = 'emc_cds.fasta'

polycistron_sequences = read_fasta(polycistron_path)
emc_sequences = read_fasta(emc_path)

polycistron_codons, total_polycistron_codons = count_codons(polycistron_sequences)
emc_codons, total_emc_codons = count_codons(emc_sequences)

polycistron_codons_percentage = calculate_percentage(polycistron_codons, total_polycistron_codons)
emc_codons_percentage = calculate_percentage(emc_codons, total_emc_codons)

df_polycistron = pd.DataFrame(list(polycistron_codons_percentage.items()), columns=['Codon', 'Group 1 Percentage'])
df_emc = pd.DataFrame(list(emc_codons_percentage.items()), columns=['Codon', 'Group 2 Percentage'])

df_combined = pd.merge(df_polycistron, df_emc, on='Codon', how='outer').fillna(0)
df_combined = df_combined.sort_values(by='Group 1 Percentage', ascending=False)  # Sort based on Group 1 data

df_combined_full_corrected_all_rows = df_combined.copy()
df_combined_full_corrected_all_rows.reset_index(drop=True, inplace=True)

df_combined_full_corrected_all_rows.head(len(df_combined_full_corrected_all_rows))
