#!/bin/bash
set -e

PREFIX="sample"

GFF_INPUT="input_annotations.gff"
FASTA_INPUT="../reference_genome.fa"
TYPES=("upstream" "downstream")  # Add any types as needed

GFF_WITH_TSS_TES="${PREFIX}_with_tss_tes.gff"
MRNA_MATCH_FILE="${PREFIX}_mRNA_gene_match.txt"
START_GFF="${PREFIX}_start_codon.gff"
STOP_GFF="${PREFIX}_stop_codon.gff"
START_BED="${PREFIX}_start_codon.bed"
STOP_BED="${PREFIX}_stop_codon.bed"
STOP_FASTA="${PREFIX}_stop.fa"

agat_sp_add_start_and_stop.pl -gff "${GFF_INPUT}" -fa "${FASTA_INPUT}" -o "${GFF_WITH_TSS_TES}"

awk '$3=="mRNA"' "${GFF_WITH_TSS_TES}" | cut -f 9 | awk -F ';' '{print $1, $2, $3}' OFS='\t' | sed 's/ID=//g; s/Name=//g; s/Parent=//g; s/-RA-00001//g' > "${MRNA_MATCH_FILE}"

awk '$3=="start_codon"' "${GFF_WITH_TSS_TES}" > "${START_GFF}"
agat_convert_sp_gff2bed.pl --gff "${START_GFF}" -o "${START_BED}"

awk '$3=="stop_codon"' "${GFF_WITH_TSS_TES}" > "${STOP_GFF}"
agat_convert_sp_gff2bed.pl --gff "${STOP_GFF}" -o "${STOP_BED}"

bedtools getfasta -fi "${FASTA_INPUT}" -fo "${STOP_FASTA}" -bed "${STOP_BED}" -s -nameOnly

for TYPE in "${TYPES[@]}"; do
    NAMES_FILE="${TYPE}_names.txt"
    TRANSCRIPT_IDS="${PREFIX}_${TYPE}_transcript_ID.txt"
    START_TYPE_BED="${PREFIX}_${TYPE}_start.bed"
    START_CLEAN_BED="${PREFIX}_${TYPE}_start_clean.bed"
    START_CLEAN_FA="${PREFIX}_${TYPE}_start_clean.fa"

    grep -w -F -f "${NAMES_FILE}" "${MRNA_MATCH_FILE}" | cut -f 1 > "${TRANSCRIPT_IDS}"
    grep -w -F -f "${TRANSCRIPT_IDS}" "${START_BED}" > "${START_TYPE_BED}"
    sort -k1,1 -k2,2n "${START_TYPE_BED}" | awk -v type="${TYPE^}" '{
        if($6=="+") {
            print $1, ($2-17), $3, type "_TSS_" NR, $5, $6
        } else if($6=="-") {
            print $1, $2, ($3+17), type "_TSS_" NR, $5, $6
        }
    }' OFS='\t' > "${START_CLEAN_BED}"
    bedtools getfasta -fi "${FASTA_INPUT}" -fo "${START_CLEAN_FA}" -bed "${START_CLEAN_BED}" -s -nameOnly
done

echo "All steps completed successfully."