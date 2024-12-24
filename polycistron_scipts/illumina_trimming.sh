# Short-reads RNAseq / ATACseq trimming
for read1 in $(find "$SOURCE" -name '*_R1_001.fastq');
 do
	read2=$(echo $read1 | sed 's/_R1_001.fastq/_R2_001.fastq/g')
	trim_galore --cores 8 --paired --fastqc "$read1" "$read2" -o "$TRIMMED"
done
for files in $(find "$SOURCE" -name '*\.fastq');
 do
	trim_galore --cores 8 \
	--fastqc -o trimmed \
	"$files"
done