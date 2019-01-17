#!bin/bash
BLAST_DATABASE=/projects/redser3-notbackedup/ftp/NCBI/BLAST/blastdb/nt
for SAMPLE in $@; do
	# Compare reads to BLAST database
	blastn -num_threads 15 -query ${SAMPLE} -db ${BLAST_DATABASE} -outfmt 6 -out ${SAMPLE}.BLAST.txt
done