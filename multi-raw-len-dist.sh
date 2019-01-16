#!bin/bash

echo "Running..."
for filename in $@; do
	echo $filename
	#Identify and unzip file
	gunzip ${filename}
	echo "Unzipped file"
	prefix=$(echo $filename | awk -F'[_]' '{print $1}')
	echo Prefix: $prefix
	outfile=$prefix-rawlength.txt
	#Iterate over sequence length before adapter sequence (R2)
	LEN=0
	while [ "$LEN" -lt 71 ]; do
		#echo Length: $LEN
		VAR=$(grep "^[A-Z]\{$LEN\}AGATCGG" ${prefix}*R2_001.fastq | wc -l)
		#echo Occurences: $VAR
		#Write values to CSV
		echo -n "$LEN	$VAR" >> $outfile
		echo "" >> $outfile
		let LEN=LEN+1
	done

	#Close file
	gzip ${prefix}*R2_001.fastq
done
