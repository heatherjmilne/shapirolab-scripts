#!bin/bash

#Identify and unzip file
### EDIT FILENAME ###
filename=HM198-1

gunzip ${filename}*R2_001.fastq.gz
outfile=$filename-rawlength.txt
#Iterate over sequence length before adapter sequence (R2)
LEN=0
while [ "$LEN" -lt 71 ]; do
	#echo Length: $LEN
	VAR=$(grep "^[A-Z]\{$LEN\}AGATCGG" ${filename}*R2_001.fastq | wc -l)
	#echo Occurences: $VAR
	#Write values to CSV
	echo -n "$LEN	$VAR" >> $outfile
	echo "" >> $outfile
	let LEN=LEN+1
done

#Close file
gzip ${filename}*R2_001.fastq
