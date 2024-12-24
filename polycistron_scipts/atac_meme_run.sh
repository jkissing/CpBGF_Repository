# ATAC-seq 

bowtie2-build <INPUT_REF.fasta> <BOWTIE2_INDEX_PREFIX>

bowtie2  --very-sensitive  -k 10  -x <BOWTIE2_INDEX_PREFIX>  \
-1 <TRIMMED_READ_1>  \
-2 <TRIMMED_READ_2>  \
-p 8 | samtools sort -@ 8 -O bam -o <ATAC_BOWTIE2.bam>

java -jar picard.jar MarkDuplicates \
-I <ATAC_BOWTIE2.bam> \
-O <ATAC_BOWTIE2_DEDUP.bam> \
-M <ATAC_BOWTIE2_DEDUP_STATS.txt>

bamPEFragmentSize -o atac_frag_length.svg \
--maxFragmentLength 1000 \
-b <ATAC_BOWTIE2_DEDUP.bam> \
-p 8 --plotFileFormat svg

java -jar picard.jar CollectInsertSizeMetrics \
-I <ATAC_BOWTIE2_DEDUP.bam> \
-O <ATAC_BOWTIE2_INSERT_STAT.txt> \
-H <ATAC_BOWTIE2_INSERT_STAT.pdf> \
-M 0.5

df <- read.table("<ATAC_BOWTIE2_INSERT_STAT.txt>", header = T, sep = "\t")

ggplot(df, aes(x=insert_size)) +
 geom_point(aes(y=Experimental), colour = "#7851A9", size = 2) +
 geom_line(aes(y=Experimental), colour = "#7851A9", alpha=0.9, size = 2) +
 geom_point(aes(y=Control), colour = "#A95156", size = 2) +
 geom_line(aes(y=Control), colour = "#A95156", alpha=0.5, size = 2) +
 scale_x_continuous(limits=c(0,1500),breaks = seq(0, 1500, by = 100)) + 
 ylim(0,200000) + ggtitle("ATAC_BGF") + theme_bw()

macs3 hmmratac \
-b <ATAC_BOWTIE2_DEDUP.bam> \
--outdir macs3_hmmr \
-n macs3_hmmr

annotatePeaks.pl <MACS3_HMMR.bed> <INPUT_REF.fasta> \
-annStats macs3_peak_annotation_stats.txt \
-gtf <INPUT_REF.gtf> > <HOMNER_OUTPUT.txt>

# MEME motif discovery

bedtools getfasta \
-fi <INPUT_REF.fasta> \
-fo <OUTPUT_REGIONS.fasta> \
-bed <MACS3_HMMR.bed> \
-s -nameOnly

meme <OUTPUT_REGIONS.fasta> \
-oc macs3_hmmr_meme \
-objfun classic \
-dna -nmotifs 3 -p 8