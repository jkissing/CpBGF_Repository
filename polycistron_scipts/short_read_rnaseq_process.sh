# Short-read RNAseq analysis

#!/bin/bash
set -e

GENOME_DIR="STAR_genome"
READS_DIR="/scratch/rx99114/cp_pt/PT_expression/0h_rnaseq"
OUTPUT_PREFIX_DIR="rnaseq/bam"
THREADS=60

mkdir -p "${OUTPUT_PREFIX_DIR}"

for read_file in "${READS_DIR}"/*_trimmed.fq; do
    base=$(basename "${read_file}" _trimmed.fq)
    output_prefix="${OUTPUT_PREFIX_DIR}/${base}_"
    
    STAR --runMode alignReads \
        --runThreadN "${THREADS}" \
        --genomeDir "${GENOME_DIR}" \
        --readFilesIn "${read_file}" \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --bamRemoveDuplicatesType UniqueIdentical \
        --twopassMode Basic \
        --outFileNamePrefix "${output_prefix}" \
        --outWigStrand Stranded \
        --outSAMstrandField intronMotif \
        --alignIntronMax 2500 \
        --chimOutType Junctions \
        --limitBAMsortRAM 1000000000
done

#!/bin/bash
set -e

INPUT_DIR="path/to/bam_files"
OUTPUT_DIR="path/to/output_bw_files"
THREADS_INDEX=20
THREADS_COVERAGE=60
STRANDS=("reverse" "forward")

mkdir -p "${OUTPUT_DIR}"

for bam in "${INPUT_DIR}"/*.bam; do
    samtools index -@ "${THREADS_INDEX}" "${bam}"
    base=$(basename "${bam}" .bam)
    for strand in "${STRANDS[@]}"; do
        suffix="REV"
        if [ "${strand}" == "forward" ]; then
            suffix="FWD"
        fi
        bamCoverage -b "${bam}" -o "${OUTPUT_DIR}/${base}_${suffix}.bw" -bs 1 -p "${THREADS_COVERAGE}" --filterRNAstrand "${strand}"
    done
done

echo "All steps completed successfully."