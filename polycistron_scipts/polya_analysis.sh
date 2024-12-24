#!/bin/bash
set -e

RAW_DIR="raw"
FQ_DIR="fq"
SEQUENCING_SUMMARY="${FQ_DIR}/sequencing_summary.txt"
ALL_FASTQ="${FQ_DIR}/all.fastq"
PREFIX="sample"
THREADS=8
BAM_FILE="${PREFIX}_aligned_sorted_hq.bam"
GENOME_FA="${PREFIX}.fa"


nanopolish index \
    --directory="${RAW_DIR}" \
    --sequencing-summary="${SEQUENCING_SUMMARY}" \
    "${ALL_FASTQ}"


POLYA_RESULTS="${PREFIX}_polya_results.tsv"

nanopolish polya \
    --threads=${THREADS} \
    --reads="${FQ_DIR}/pass.fastq" \
    --bam="${BAM_FILE}" \
    --genome="${GENOME_FA}" \
    > "${POLYA_RESULTS}"


PASS_READNAMES="${PREFIX}_pass_readname.txt"
PASS_POLYA_LENGTH="${PREFIX}_pass_polya_length.txt"

grep 'PASS' "${POLYA_RESULTS}" | cut -f1 | sort -n > "${PASS_READNAMES}"
grep 'PASS' "${POLYA_RESULTS}" | cut -f1,9 | sort -k1n | cut -f2 > "${PASS_POLYA_LENGTH}"


PASS_FASTQ="${FQ_DIR}/pass.fastq"
PASS_FASTA="../${PREFIX}_pass.fasta"
PASS_READ_FASTA="${PREFIX}_pass_read.fasta"

seqkit fq2fa "${PASS_FASTQ}" -o "${PASS_FASTA}"
seqkit grep -f "${PASS_READNAMES}" "${PASS_FASTA}" -o "${PASS_READ_FASTA}"

samtools faidx "${PASS_READ_FASTA}"

PASS_READ_FASTA_FAI="${PASS_READ_FASTA}.fai"
PASS_POLYA_TRANSCRIPT_LENGTH="${PREFIX}_pass_polya_transcript_length.txt"

sort -k1n "${PASS_READ_FASTA_FAI}" | cut -f2 > "${PASS_POLYA_TRANSCRIPT_LENGTH}"

TRANSCRIPT_VS_POLYA="${PREFIX}_transcript_vs_polya_length.txt"
FINAL_OUTPUT="${PREFIX}_transcript_VS_polya_length.txt"

paste "${PASS_POLYA_LENGTH}" "${PASS_POLYA_TRANSCRIPT_LENGTH}" > "${TRANSCRIPT_VS_POLYA}"
awk '{printf "%.0f\t%.0f\n", $1, $2}' "${TRANSCRIPT_VS_POLYA}" > "${FINAL_OUTPUT}"
