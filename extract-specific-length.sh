#!bin/bash

size=42
echo "Running..."
# 1. Extract all of the reads of length 42 prior to the adapter sequence
# 2. Convert to FASTA

echo "Extracting reads of length $size..."
for filename in $@; do
	echo $filename
	#Identify and unzip file
	gunzip ${filename}
	echo "Unzipped file"
	prefix=$(echo $filename | awk -F'[_]' '{print $1}')
	echo Prefix: $prefix
	outfile=${prefix}-${size}bp-R2_001.fastq
	#Search for all lines that are a particular length prior to the adapter sequence
	#Output to the outfile in FASTQ format
	grep -B 1 -A 2 "^[A-Z]\{$size\}AGATCGG" ${prefix}*.fastq > ${outfile}
	sed -i '' '/--/d' ${outfile}
	gzip ${prefix}*_R2_001.fastq
	#Now the file is in FASTQ - convert to FASTA
	/soe/pheintzman/bin/fastx_toolkit-0.0.13.2/src/fastq_to_fasta/fastq_to_fasta -n -Q33 -r -i ${outfile} -o ${prefix}-${size}-42bp-R2_001.fasta
	#gzip ${outfile}
	#Close file
done

echo "Creating multiple sequence alignment"


# 3. Multiple sequence alignment / Motif finder
# 4. BLAST the NR database
# If the results are weird or bad because of the adapter sequences, then trim them and try again.
