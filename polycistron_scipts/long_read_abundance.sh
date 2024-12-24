#!/bin/bash
set -e

PREFIX="cp"
INPUT_PROT="cp_pre_proteomic.txt"
TRANSCRIPTOME_FA="cp_t2tbgfname_transcriptome_all_wIUTR.fa"
REFERENCE_FA="../cp_bgf.fa"
GENOME_DB="bgf"
GENOME_GENOME="../cp_bgf.genome"
OUTPUT_DIR="blast_results"
FASTA_DIR="fasta_files"
MAP_DIR="mapping_results"
THREADS_MINIMAP2=8
THREADS_BLAST=12

mkdir -p "${OUTPUT_DIR}" "${FASTA_DIR}" "${MAP_DIR}"

declare -A ID_TYPES=(
    ["IUM"]="pt_all_IUM_id.txt"
    ["IDM"]="pt_all_IDM_id.txt"
)

for TYPE in "${!ID_TYPES[@]}"; do
    ID_FILE="${ID_TYPES[$TYPE]}"
    OUTPUT_FA="${FASTA_DIR}/${PREFIX}_pt_all_${TYPE}.fa"
    seqkit grep -n -f "${ID_FILE}" "${TRANSCRIPTOME_FA}" > "${OUTPUT_FA}"
done

ALL_FA="${FASTA_DIR}/${PREFIX}_prot_all.fa"
tail -n+2 "${INPUT_PROT}" | cut -f5 | awk '{i++; printf(">Seq%d\n%s\n", i, $1)}' > "${ALL_FA}"

makeblastdb -dbtype nucl -in "${REFERENCE_FA}" -title "${GENOME_DB}" -out "${GENOME_DB}"

declare -A MAPPING_CONFIGS=(
    ["isoseq_IUM"]="minimap2 -x splice:hq --for-only -G 2500"
    ["isoseq_IDM"]="minimap2 -x splice:hq --for-only -G 2500"
    ["drs_arizona_IUM"]="minimap2 -x splice -uf -k14 --for-only -G 2500"
    ["drs_arizona_IDM"]="minimap2 -x splice -uf -k14 --for-only -G 2500"
    ["drs_inrae_IUM"]="minimap2 -x splice -uf -k14 --for-only -G 2500"
    ["drs_inrae_IDM"]="minimap2 -x splice -uf -k14 --for-only -G 2500"
)

declare -A READ_FILES=(
    ["isoseq_IUM"]="/scratch/rx99114/isoseq_bgf/isoseq_cp.fq"
    ["isoseq_IDM"]="/scratch/rx99114/isoseq_bgf/isoseq_cp.fq"
    ["drs_arizona_IUM"]="../../drs_arizona/arizona_drs_dorado.fq"
    ["drs_arizona_IDM"]="../../drs_arizona/arizona_drs_dorado.fq"
    ["drs_inrae_IUM"]="../drs_control_treated_combined.fq"
    ["drs_inrae_IDM"]="../drs_control_treated_combined.fq"
)

for CONFIG in "${!MAPPING_CONFIGS[@]}"; do
    TOOL_CMD="${MAPPING_CONFIGS[$CONFIG]}"
    READ_FILE="${READ_FILES[$CONFIG]}"
    QUERY_FA="${FASTA_DIR}/${PREFIX}_pt_all_${CONFIG##*_}.fa"
    OUTPUT_PAF="${MAP_DIR}/${PREFIX}_${CONFIG}_map.paf"
    ${TOOL_CMD} -t "${THREADS_MINIMAP2}" "${QUERY_FA}" "${READ_FILE}" > "${OUTPUT_PAF}"
done

declare -A BLAST_QUERIES=(
    ["sporozoite"]="fasta_files/cp_pt_all_IUM.fa 1e-10 5"
    ["all"]="fasta_files/cp_prot_all.fa 1e-4 10"
)

for TYPE in "${!BLAST_QUERIES[@]}"; do
    IFS=' ' read -r QUERY_FA EVALUE ALIGNMENTS <<< "${BLAST_QUERIES[$TYPE]}"
    OUT_TXT="${OUTPUT_DIR}/${PREFIX}_bgf_prot_aligned_${TYPE}.txt"
    tblastn -query "${QUERY_FA}" -db "${GENOME_DB}" \
        -outfmt 6 -max_intron_length 28000 \
        -out "${OUT_TXT}" -evalue "${EVALUE}" \
        -num_alignments "${ALIGNMENTS}" -num_threads "${THREADS_BLAST}"
done

cat "${OUTPUT_DIR}/${PREFIX}_bgf_prot_aligned_sporozoite.txt" "${OUTPUT_DIR}/${PREFIX}_bgf_prot_aligned_all.txt" > "${OUTPUT_DIR}/${PREFIX}_bgf_prot_aligned.txt"

awk '{if($9<$10) print $2, "tblastn", "match", $9, $10, ".", "+", ".", "ID=Peptide_"$1; else if($9>$10) print $2, "tblastn", "match", $10, $9, ".", "-", ".", "ID=Peptide_"$1}' OFS='\t' "${OUTPUT_DIR}/${PREFIX}_bgf_prot_aligned.txt" | sort -k1,1 -k4,4n > "${OUTPUT_DIR}/${PREFIX}_bgf_prot_aligned.gff"

awk '{if($9<$10) print $2, "tblastn", "match", $9, $10, $11, "+", ".", "ID=Peptide_"$1";Target="$1" "$7" "$8; else if($9>$10) print $2, "tblastn", "match", $10, $9, $11, "-", ".", "ID=Peptide_"$1";Target="$1" "$7" "$8}' OFS='\t' "${OUTPUT_DIR}/${PREFIX}_bgf_prot_aligned.txt" | sort -k1,1 -k4,4n > "${OUTPUT_DIR}/${PREFIX}_bgf_prot_aligned_full.gff"

awk '{if($9<$10) print $2, $9, $10, $1, $12, "+"; else if($9>$10) print $2, $10, $9, $1, $12, "-"}' OFS='\t' "${OUTPUT_DIR}/${PREFIX}_bgf_prot_aligned.txt" | sort -k1,1 -k2,2n > "${OUTPUT_DIR}/${PREFIX}_bgf_prot_aligned.bed"

bedToBam -i "${OUTPUT_DIR}/${PREFIX}_bgf_prot_aligned.bed" -g "${GENOME_GENOME}" > "${OUTPUT_DIR}/${PREFIX}_bgf_prot_aligned.bam"


# Counting PT read supports
declare -A COUNT_FILES=(
    ["isoseq_IUM"]="mapping_results/isoseq_pt_all_IUM_map.paf"
    ["isoseq_IDM"]="mapping_results/isoseq_pt_all_IDM_map.paf"
    ["drs_arizona_IUM"]="mapping_results/drs_arizona_pt_all_IUM_map.paf"
    ["drs_arizona_IDM"]="mapping_results/drs_arizona_pt_all_IDM_map.paf"
    ["drs_inrae_IUM"]="mapping_results/drs_inrae_pt_all_IUM_map.paf"
    ["drs_inrae_IDM"]="mapping_results/drs_inrae_pt_all_IDM_map.paf"
)
COUNT_OUTPUT="${OUTPUT_DIR}/pt_counts.tsv"
echo -e "Source\tType\tCount" > "${COUNT_OUTPUT}"

for KEY in "${!COUNT_FILES[@]}"; do
    FILE="${COUNT_FILES[$KEY]}"
    if [[ "${KEY}" == isoseq* ]]; then
        CONDITION="Isoseq"
    else
        CONDITION="DRS_${KEY%%_*}"
    fi
    TYPE="${KEY##*_}"
    if [[ "${TYPE}" == "IUM" || "${TYPE}" == "IDM" ]]; then
        if [[ "${KEY}" == isoseq_* ]]; then
            FILTER=">1.2"
        else
            FILTER=">1.2"
        fi
    fi
    COUNT=$(cut -f6,2,7 "${FILE}" | awk '{print $2,($1/$3)}' OFS='\t' | awk '$2>1.2' | cut -f1 | sort -V | uniq -c | awk '{print $2"\t"$1}')
    echo -e "${CONDITION}\t${TYPE}\t${COUNT}" >> "${COUNT_OUTPUT}"
done

for KEY in "${!COUNT_FILES[@]}"; do
    FILE="${COUNT_FILES[$KEY]}"
    if [[ "${KEY}" == isoseq_* ]]; then
        CONDITION="Isoseq"
    else
        CONDITION="${KEY%%_*}_DRS"
    fi
    TYPE="${KEY##*_}"
    COUNT=$(cut -f6,2,7 "${FILE}" | awk '{print $2,($1/$3)}' OFS='\t' | awk '$2<=1' | cut -f1 | sort -V | uniq -c | awk '{print $2"\t"$1}')
    echo -e "${CONDITION}\t${TYPE}\t${COUNT}" >> "${COUNT_OUTPUT}"
done
