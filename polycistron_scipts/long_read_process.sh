# Long-read RNA-seq processing and polycistron annotation
lima <RAW_CCS.bam> \
primers.fasta <LIMA.bam> \
--isoseq --peek-guess \
--dump-removed \
-j 64
isoseq3 refine \
<LIMA.bam> \
primers.fasta <FLNC.bam> \
-j 64 --require-polya -v
pbmm2 index <INPUT_REF.fasta> <OUTPUT_REF.mmi>
pbmm2 align \
<OUTPUT_REF.mmi> <FLNC.bam> <ALIGNED.bam> \
--sort -j 32 -J 32 --preset ISOSEQ \
-G 2500 --log-level INFO
dorado basecaller rna002_70bps_hac@v3 <POD5_INPUT_FOLDER> \
-v --emit-fastq  \
> <OUTPUT_DRS.fastq>

minimap2 -ax splice -uf -k14 -t 8 \
-G 2500 \
<REF_GENOME.fasta> <DRS.fastq> > <DRS.sam>

samtools view -h -@ 8 -S -b <DRS.sam> \
| samtools sort -@ 8 -o <DRS_SORTED.bam>
samtools index -@ 8  <DRS_SORTED.bam>

