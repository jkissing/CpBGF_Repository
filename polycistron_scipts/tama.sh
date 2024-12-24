# TAMA Collapse
tama_collapse.py -s <GENOME_ALIGNED_LONGREAD_RNASEQ.sam> \
-f <INPUT_REF.fasta> \
-p tama_high_quality_collapse \
-x no_cap \
-d merge_dup \
-a 100 -z 100 \
-sj sj_priority -lde 1 -sjt 20


# TAMA Refine
tama_merge.py -f tama_merge_filelist.txt -p tama_merge -d merge_dup
python ~/tama_go/read_support/tama_read_support_levels.py \
-f filelist.txt -o tama_merge -m tama_merge_merge.txt
bedtools getfasta -name -split -s \
-fi <reference genome in fa> \
-bed <tama collapse in bed> \
-fo <output in fa format>
python tama_go/orf_nmd_predictions/tama_orf_seeker.py \
-f <output in fa format> -o <orf output in fa format>
blastp -evalue 1e-10 \
-num_threads 20 \
-db /db/uniprot/latest/uniref100 \
-query <orf output in fa format> \
-out <blastp output in txt format>
python tama_go/orf_nmd_predictions/tama_orf_blastp_parser.py \
-b <blastp output in txt format> -o <blastp output intermediate in txt format>
python tama_go/orf_nmd_predictions/tama_cds_regions_bed_add.py \
-p <blastp output intermediate in txt format> \
-a <tama collapse in bed> \
-f <collapse sequence output in fa format> \
-o <final CDS added output in bed format>
-f <final CDS added output in bed format> -o <fragment removed in bed format, without .bedt> -cds longest_cds \
-id original_id
python tama_go/format_converter/tama_format_id_filter.py \
-b <final output in bed format> -o <final sorted output in bed format> \
-s custom -r 2,1
gt bed_to_gff3 -thicktype CDS \
-featuretype mRNA -blocktype exon \
<final sorted output in bed format> > <final sorted output in gff format>
python tama_go/format_converter/tama_convert_bed_gtf_ensembl_orf_nmd.py \
<fragment removed in bed format> <fragment removed in gtf format>
agat_convert_sp_gxf2gxf.pl -g <fragment removed in gtf format> \
-o <final in gff format>