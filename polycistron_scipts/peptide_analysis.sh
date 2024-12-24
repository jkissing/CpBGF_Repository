#!/bin/bash
set -e

PREFIX="cp"
INPUT_PROT="cp_pre_proteomic.txt"
REFERENCE_FA="../cp_bgf.fa"
GENOME_DB="bgf"
GENOME_GENOME="../cp_bgf.genome"
OUTPUT_DIR="blast_results"
FASTA_DIR="fasta_files"

mkdir -p "${OUTPUT_DIR}" "${FASTA_DIR}"

declare -A FILTERS=(
    ["sporozoite"]="Sporozoite"
    ["other"]="!Sporozoite"
    ["all"]=""
)

for TYPE in "${!FILTERS[@]}"; do
    FILTER=${FILTERS[$TYPE]}
    if [ "$TYPE" == "all" ]; then
        awk 'NR>1' "${INPUT_PROT}" | cut -f5 | awk '{i++; printf(">Seq%d\n%s\n", i, $1)}' > "${FASTA_DIR}/${PREFIX}_${TYPE}_raw.fa"
    else
        awk -F '\t' -v filter="$FILTER" '$3 ~ filter' "${INPUT_PROT}" | cut -f5 | awk '{i++; printf(">Seq%d\n%s\n", i, $1)}' > "${FASTA_DIR}/${PREFIX}_${TYPE}_raw.fa"
    fi
done

makeblastdb -dbtype nucl -in "${REFERENCE_FA}" -title "${GENOME_DB}" -out "${GENOME_DB}"

declare -A BLAST_PARAMS=(
    ["sporozoite"]="1e-10 5"
    ["all"]="1e-4 10"
)

for TYPE in "sporozoite" "all"; do
    QUERY_FA="${FASTA_DIR}/${PREFIX}_${TYPE}_raw.fa"
    OUT_TXT="${OUTPUT_DIR}/${PREFIX}_bgf_prot_aligned_${TYPE}.txt"
    EVALUE=$(echo "${BLAST_PARAMS[$TYPE]}" | awk '{print $1}')
    ALIGNMENTS=$(echo "${BLAST_PARAMS[$TYPE]}" | awk '{print $2}')
    tblastn -query "${QUERY_FA}" -db "${GENOME_DB}" \
        -outfmt 6 -max_intron_length 28000 \
        -out "${OUT_TXT}" -evalue "${EVALUE}" \
        -num_alignments "${ALIGNMENTS}" -num_threads 12
done

cat "${OUTPUT_DIR}/${PREFIX}_bgf_prot_aligned_sporozoite.txt" "${OUTPUT_DIR}/${PREFIX}_bgf_prot_aligned_all.txt" > "${OUTPUT_DIR}/${PREFIX}_bgf_prot_aligned.txt"
awk '{if($9<$10) print $2, "tblastn", "match", $9, $10, ".", "+", ".", "ID=Peptide_"$1; else if($9>$10) print $2, "tblastn", "match", $10, $9, ".", "-", ".", "ID=Peptide_"$1}' OFS='\t' "${OUTPUT_DIR}/${PREFIX}_bgf_prot_aligned.txt" | sort -k1,1 -k4,4n > "${OUTPUT_DIR}/${PREFIX}_bgf_prot_aligned.gff"
awk '{if($9<$10) print $2, "tblastn", "match", $9, $10, $11, "+", ".", "ID=Peptide_"$1";Target="$1" "$7" "$8; else if($9>$10) print $2, "tblastn", "match", $10, $9, $11, "-", ".", "ID=Peptide_"$1";Target="$1" "$7" "$8}' OFS='\t' "${OUTPUT_DIR}/${PREFIX}_bgf_prot_aligned.txt" | sort -k1,1 -k4,4n > "${OUTPUT_DIR}/${PREFIX}_bgf_prot_aligned_full.gff"
awk '{if($9<$10) print $2, $9, $10, $1, $12, "+"; else if($9>$10) print $2, $10, $9, $1, $12, "-"}' OFS='\t' "${OUTPUT_DIR}/${PREFIX}_bgf_prot_aligned.txt" | sort -k1,1 -k2,2n > "${OUTPUT_DIR}/${PREFIX}_bgf_prot_aligned.bed"
bedToBam -i "${OUTPUT_DIR}/${PREFIX}_bgf_prot_aligned.bed" -g "${GENOME_GENOME}" > "${OUTPUT_DIR}/${PREFIX}_bgf_prot_aligned.bam"
