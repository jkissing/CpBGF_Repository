# Comands Used:

•	ONT gDNA basecalling:

o	Dorado v0.5.0:
dorado basecaller \
-r --emit-fastq \
dna_r9.4.1_e8_sup@v3.6 \
<input_fast5> > CpBGF_ONT_raw.fastq

•	Chimeric read detection:

o	Minimap2/yacrd:
minimap2 -x ava-ont -g 500 CpBGF_ONT_raw.fastq CpBGF_ONT_raw.fastq > overlap.paf
yacrd -i overlap.paf -o reads.yacrd

•	Genome Assembly:

o	NECAT/0.0.1:
necat.pl correct necat_config.txt
necat.pl assemble necat_config.txt
necat.pl bridge necat_config.txt

o	necat_config.txt:
PROJECT=BGF_NECAT
ONT_READ_LIST=BGF_ONT.txt
GENOME_SIZE=9200000
THREADS=5
MIN_READ_LENGTH=3000
PREP_OUTPUT_COVERAGE=30
OVLP_FAST_OPTIONS=-n 500 -z 20 -b 2000 -e 0.5 -j 0 -u 1 -a 1000
OVLP_SENSITIVE_OPTIONS=-n 500 -z 10 -e 0.5 -j 0 -u 1 -a 1000
CNS_FAST_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0
CNS_SENSITIVE_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0
TRIM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 1 -a 400
ASM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400
NUM_ITER=2
CNS_OUTPUT_COVERAGE=30
CLEANUP=1
USE_GRID=false
GRID_NODE=0
SMALL_MEMORY=0
POLISH_CONTIGS=true

•	Genome Polishing:

o	NextPolish v.1.2.4 :
sgs_options = -max_depth 200, lgs_options = -min_read_len 1k -max_read_len 100k, and lgs_minimap2_options = -x map-ont.

•	Genome Annotation:

o	BRAKER2:
		braker.pl --species=yourSpecies --genome=genome.fasta \
       		--rnaseq_sets_ids=SRA_ID1,SRA_ID2 \
       		--rnaseq_sets_dirs=/path/to/local/fastq/files/ \
		--UTR=on --addUTR=on --ab_initio

o	Liftoff:
liftoff -g CpIOWA-ATCC_reference.gff -o CpBGF_liftoff.gff -copies -cds CpBGF_genome.fasta CpIOWA-ATCC_genome.fasta


•	Iso-Seq analysis:

o	SMRTlink v10.1:
isoseq3 refine \
Isoseq3_ccs.fl.NEB_5p--NEB_Clontech_3p.bam \
primers.fasta cp.flnc.bam \
-j 8 --require-polya -v
 		pbmm2 index cp_bgf.fa cp_bgf.mmi
pbmm2 align \
cp_bgf.mmi cp.flnc.bam cp_bgf_isoseq_aligned.bam \
--sort -j 8 -J 8 --preset ISOSEQ \
-G 2500 --log-level INFO

o	TAMA:
tama_collapse.py -s tama/ cp_bgf_isoseq_aligned.sam \
-f cp_bgf.fa -p tama_cp_bgf \
-x no_cap -d merge_dup \
-a 100 -z 100 \
-sj sj_priority -lde 1 -sjt 20

•	BUSCO run:
busco -i <genome.fasta> -l ./apicomplexa_odb10 -m genome -o <output>

•	Short-read RNAseq mapping:

o	STAR v2.7.10:
STAR --runMode alignReads \
--runThreadN 16 \
--genomeDir STAR_genome \
--readFilesIn <input_read_1.fq> <input_read_2.fq> \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--bamRemoveDuplicatesType UniqueIdentical \
--twopassMode Basic \
--outFileNamePrefix <output_name> \
--outWigStrand Stranded \
--outSAMstrandField intronMotif \
--alignIntronMax 2500 \
--chimOutType Junctions \
--limitBAMsortRAM 1000000000

•	Protein sequence extraction:

o	AGAT v0.9.2:
agat_sp_extract_sequences.pl -p \
-g <input.gff> -f <genome.fa> \
-o <protein_output.fa>

•	Starting amino acid analysis:

o	Seqkit v2.5.1:
seqkit fx2tab <protein_input.fa> \
| cut -f 2 | awk '{print substr($1,1,1)}' \
| sort -V | uniq -c > <output.txt>

•	Variant Call:

o	BWA/SAMTOOLS:
bwa mem -v 2 -M -t 10 $GENOME $R1 $R2 > $SAMPLE.sam
samtools view -b -S -o $SAMPLE.bam $SAMPLE.sam

o	PICARD:
java -Xmx2g -classpath "$PICARD" -jar $SortSam INPUT=$SAMPLE.bam OUTPUT=$SAMPLE.s.bam VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate
java -Xmx2g -classpath "$PICARD" -jar $MarkDuplicates INPUT=$SAMPLE.s.bam OUTPUT=$SAMPLE.sd.bam METRICS_FILE=$SAMPLE.dedup.metrics REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true
java -Xmx2g -classpath "$PICARD" -jar $ARreadgroups INPUT=$SAMPLE.sd.bam OUTPUT=$SAMPLE.sdr.bam SORT_ORDER=coordinate RGID=$SAMPLE RGLB=$SAMPLE RGPL=illumina RGPU=$SAMPLE RGSM=$SAMPLE VALIDATION_STRINGENCY=LENIENT
java -Xmx2g -classpath "$PICARD" -jar $BuildBamIndex INPUT=$SAMPLE.sdr.bam VALIDATION_STRINGENCY=LENIENT

o	GATK:
gatk HaplotypeCaller -R $GENOME -I $SAMPLE.sdrsm.bam -O $SAMPLE.GATK.vcf -ploidy $PLOIDY -stand-call-conf 30
gatk VariantFiltration -R $GENOME -V $SAMPLE.GATK.vcf -O $SAMPLE.GATK_Filtered.vcf -cluster 3 -window 10 -filter "QUAL < 30.0 || DP < 10 || QD < 1.5 || MQ < 25.0" --filter-name "StdFilter" -filter "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filter-name "HARD_TO_VALIDATE" -filter "MQ < 40.0 || FS > 60.0" --filter-name "gatkFilter"

Other commands and pipelines used in data processing were executed using their corresponding default parameters. Bedtools is used for calculating the number of basepairs contributing to annotations through merging genomic regions with overlapping gene features.


