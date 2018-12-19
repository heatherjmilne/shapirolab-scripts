#!/bin/bash
#######################################################################################
#	"MONICA": ILLUMMINA ANCIENT DNA ANALYSIS PIPELINE 
#######################################################################################
#	created by Dr Pete Heintzman
#	optimized by (soon to be Dr) Josh Kapp and Alisa Vershinina
# 	Last time updated: 8 June 2018 
#	Contact: vershininaalissa@gmail.com
#######################################################################################

#######################################################################################
#	DISCLAMER
#	This pipeline is written in bash.
# 	It is a collection of "for" loops that operate on different programs and scripts. 
#	The authors are not responsible for the consequences of use of this pipeline. 
#	It is a responsibility of a user to understand how this pipeline works.
#	This script is provided "as is", so use it at your own risk. 
#	User is free to change anything in this version. 
#	This version processes only PE short Illumina reads. 
#	It is named "Monica" because bioinformatic tools lack female names.
#	And we already have Jared, so I figured Monica would be a good addition.
#######################################################################################

#######################################################################################
#   This version of the pipeline has been MODIFIED from the original.
#   This is HEATHER V3 7/2/18
#   Modifications were made by Heather Milne on 7/5/2018.
#   This version is designed specifically to be used with Rachel Turba's clam samples,
#   for which I have 3 nuclear genomes each and 3 mitogenomes. It doesn't make sense
#   to run the entire script 3 times.
#   This is the MITO MAPPING script, which runs MIA.
#   Changes to this version:
#   1. Each program (ie MapDamage, blastn) run outputs its own log file for easier error handling.
#   2. blastn outputs in tabular format (-outfmt 6) for use with MEGAN software
#   3. Corrected error in the check of number of arguments, and moved it to the beginning of the file
#      before assigning variable names to arguments
#   4. Tests for the existence of each input directory before beginning. This will identify
#      typos right away, instead of waiting for the whole script to run first.

#######################################################################################
#	USAGE
#	This version runs only on edser2. Are you on edser2?
#	Before running, make this script executable:
#	chmod +x monica.sh
#	To run the pipeline:
#	./monica.sh -p=[prefix] -dir=[fastq_path] -t=[taxon] -ref=[ref_fasta_path] -refname=[ref_name] -mref=[mito_fasta_path] -mname=[mito_name] 
#		-p: extension for filenames of interest. Letters\numbers only.
#			To analyse fastq files named as AV125_S1.filename.fastq.gz, AV251_S6.filename.fastq.gz, etc - use AV.
#			To analyse only AV125_S1.filename.fastq.gz, and AV126_S7.filename.fastq.gz - use AV1 or AV12. And so on.
#			To run this pipeline on something other than our lab MiSeq data, change {SAMPLE}_L001_R1\R2_001.fastq pattern
#			to an appropriate pattern of interest (line ~187) .
#		-dir: full path to the folder with raw data (fastq.gz). Path only.
#			Example: /projects/redser2/raw/171018_M00160_0050_000000000-BFNHV (no closing slash, no filename)
#		-t: your organism you map into. Example: Bison. Can be anything.
#		-ref: full path to a reference genome in fasta format. Full path to fasta only. Reference should have bwa-compatible index. 
#			Example: /projects/redser2/genomes/Bison1.0/GCA_000754665.1_Bison_UMD1.0_genomic.fna
#		-refname: name of the reference genome. Example: UMD1.0. Can be anything.
#		-mref: full path to mito-reference genome in fasta format. Full path to fasta only.
#			Example: /projects/redser2/mito/Bison1.0/GCA_blah.fa
#		-mname: name of the reference mitochondria. Example: bison_mito. Can be anything.
#	If you analyse *capture* data, change MIA processing pipeline at the end of this script.
#######################################################################################
echo "Number of arguments: $#"
# test if we have enough arguments, exit with warning if we don't

if [ $# -lt 4 ] ### There was a TYPO here!! changed to -lt from -ge
then
	echo "Not enough arguments!"
	echo "Usage $0 -p=[prefix] -dir=[fastq_path] -t=[taxon] -ref=[ref_fasta_path] -refname=[ref_name] -mref=[mito_fasta_path] -mname=[mito_name]"
	exit 1
fi

DATE=`date +%Y-%m-%d`
echo "Today is $DATE"

# Assign variable names to arguments and tests if each directory exists
# If a directory doesn't exist, it exits the program.
for i in "$@"
do
case $i in
	-p=*|--prefix=*)
	PREFIX="${i#*=}"
	echo "PREFIX=${PREFIX}"
	shift
	;;
	# -dir=*|--raw-read-dir=*)
	# RAW_READ_DIRECTORY="${i#*=}"
	# if [ -e "$RAW_READ_DIRECTORY" ]
	# then
 #    	echo "RAW READ DIRECTORY=${RAW_READ_DIRECTORY}"
	# else
 #    	echo "RAW READ DIRECTORY $RAW_READ_DIRECTORY does not exist, exiting script"
 #    	exit 1
	# fi
	# shift
	# ;;
	-t=*|--taxon=*)
	TAXON="${i#*=}"
	echo "TAXON=$TAXON"
	shift
	;;
	# -ref=*|--ref-fasta=*)
	# REFERENCE_SEQUENCE="${i#*=}"
	# if [ -e "$REFERENCE_SEQUENCE" ]
	# then
 #    	echo "REFERENCE SEQUENCE=$REFERENCE_SEQUENCE"
	# else
 #    	echo "REFERENCE SEQUENCE $REFERENCE_SEQUENCE does not exist, exiting script"
 #    	exit 1
	# fi
	# shift
	# ;;
	# -refname=*|--ref-name=*)
	# REFERENCE_NAME="${i#*=}"
	# echo "REFERENCE NAME=$REFERENCE_NAME"
	# shift
	# ;;
	-mref=*|--mito-ref-fasta=*)
	MIA_REFERENCE_SEQUENCE="${i#*=}"
	if [ -e "$MIA_REFERENCE_SEQUENCE" ]
	then
    	echo "MIA REFERENCE SEQUENCE=$MIA_REFERENCE_SEQUENCE"
	else
    	echo "MIA REFERENCE SEQUENCES $MIA_REFERENCE_SEQUENCE does not exist, exiting script"
    	exit 1
	fi
	shift
	;;
	-mname=*|--mito-ref-name=*)
	MIA_REFERENCE_NAME="${i#*=}"
	echo "MIA REFERENCE NAME=$MIA_REFERENCE_NAME"
	shift
	;;
esac
done


# mkdir ${PWD}/${PREFIX}_shotgun_data_processing_${DATE} # output is created in the same folder, where you run the script.

# General envelopes
PROCESSING_OUTPUT=${PWD}/${PREFIX}_shotgun_data_processing_${DATE}
GET_SAMPLES=/projects/redser3-notbackedup/projects/alisa_beringia/scripts/get_samplelist.sh
CALC_STATS=/projects/redser3-notbackedup/projects/alisa_beringia/scripts/calculate_library_stats.R


### Initial setting up of directories

cd ${RAW_READ_DIRECTORY}
# bash ${GET_SAMPLES} ${PREFIX} f > ${PROCESSING_OUTPUT}/${PREFIX}-sample-list-${DATE}.txt

cd ${PROCESSING_OUTPUT}

# mkdir Sample_lists_and_progress_files
# mkdir Raw_data_symlinks
# mkdir SeqPrep_output
# mkdir BWA_analyses
# mkdir MapDamage_output
# mkdir MEGAN_analyses
mkdir MIA_analyses



# SeqPrep envelopes

SEQPREP_OUTPUT=${PROCESSING_OUTPUT}/SeqPrep_output	    # A directory to store SeqPrep output files
SEQPREP_MIN_LENGTH=30                                 	# Removes unmappably short reads. This is a fairly low cut-off
SEQPREP_OVERLAP=15                                 	    # Allows for confident merging of reads. Can be reduced to 10 if needed.
SEQPREP_LOCATION=/soe/pheintzman/bin/SeqPrep2-master	# To find SeqPrep, if using edser2

# echo "Minimum read length is set to $SEQPREP_MIN_LENGTH."
# echo "Merging reads overlapping at $SEQPREP_OVERLAP bases."
# BWA and SAMtools envelopes

BWA_OUTPUT=${PROCESSING_OUTPUT}/BWA_analyses	                                                # A directory to store BWA intermediate files and output   
INDEX_ALGORITHM=bwtsw	                                                                        # If reference is <2Gb use 'is', if >2Gb use 'bwtsw'
SEED_DISABLE=1024	                                                                            # Following ancient DNA data processing protocols
BWA_THREADS=15	                                                                                # To speed up analysis
BAM_MIN_QUALITY=20	                                                                            # Provides a fairly low cutoff

# echo "Seed disabled for seedless mapping."
# echo "Minimum mapping quality is set to $BAM_MIN_QUALITY phred."


# mapDamage envelopes - make sure REFERENCE_SEQUENCE= REFERENCE_NAME and BWA_OUTPUT are enabled

MAPDAMAGE2=/soe/pheintzman/bin/mapDamage-master/bin/mapDamage	        # To assess damage rates in the dataset
MAX_MISINCORP_FREQUENCY=0.3	                                           	# Use 0.3 if not too badly damaged, use 0.5 if badly damaged
READ_PLOT_LENGTH=25	                                               	    # The number of nucleotides to plot at the 5' and 3' ends of the read
MAX_READ_LENGTH=150	                                           	        # The maximum read length to consider

# echo "Nucleotide misincorporation frequency set to ${MAX_MISINCORP_FREQUENCY}."

# MEGAN envelopes

MEGAN_OUTPUT=${PROCESSING_OUTPUT}/MEGAN_analyses                           # A directory to store MEGAN files and output   
BLAST_DATABASE=/projects/redser3-notbackedup/ftp/NCBI/BLAST/blastdb/nt     # Location of the BLAST database 
MEGAN_THREADS=15                                                           # To speed up analysis


# MIA envelopes

MIA_OUTPUT=${PROCESSING_OUTPUT}/MIA_analyses
MIA_COVERAGE_FILTER=/projects/redser3-notbackedup/projects/pheintzman/Scripts/mia_consensus_coverage_filter.pl
FASTX_TOOLKIT=/soe/pheintzman/bin/fastx_toolkit-0.0.13.2/src
MIA_COVERAGE_FILTER_ANDRE=/projects/redser3-notbackedup/projects/common_jobs/coverage_filter_3.pl

# Prinseq envelopes

PRINSEQ_LITE=/projects/redser3-notbackedup/projects/pheintzman/Scripts/prinseq-lite.pl
PRINSEQ_GRAPHS=/projects/redser3-notbackedup/projects/pheintzman/Scripts/prinseq-graphs.pl
PRINSEQ_STATS=${PROCESSING_OUTPUT}/PRINSEQ_stats
COMPLEXITY_METHOD=dust			    # dust is the standard used by BLAST. The entropy method is an alternative, but is not widely used.
COMPLEXITY_THRESHOLD=7			    # Pretty low, but follows the PRINSEQ_LITE manual. Recommends 70 if using entropy.
COMBINE_PAIRED_END_READS=/projects/redser3-notbackedup/projects/pheintzman/Scripts/combinePairedEndReads.pl
SPLIT_PAIRED_END_READS=/projects/redser3-notbackedup/projects/pheintzman/Scripts/splitPairedEndReads.pl

# Other envelopes

GET_INSERT_SIZE=/projects/redser3-notbackedup/projects/pheintzman/Scripts/getinsertsize.py

# ### MEGAN and MIA -- Initial data processing
 
for SAMPLE in $(cat ${PREFIX}-sample-list-${DATE}.txt)
do

# Concatenate merged and unmerged reads - MEGAN and MIA
	cat ${SEQPREP_OUTPUT}/${SAMPLE}_merged.complexity_filtered.fastq ${SEQPREP_OUTPUT}/${SAMPLE}_R1_unmerged.complexity_filtered.fastq ${SEQPREP_OUTPUT}/${SAMPLE}_R2_unmerged.complexity_filtered.fastq > ${SEQPREP_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.fastq
	wait

# Convert FASTQ to FASTA - MEGAN and MIA
	${FASTX_TOOLKIT}/fastq_to_fasta/fastq_to_fasta -n -Q33 -r -i ${SEQPREP_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.fastq -o ${SEQPREP_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.fasta
	wait

# Collapse identical reads - forward, reverse, and 5' duplicates - MEGAN and MIA
	perl ${PRINSEQ_LITE} -fasta ${SEQPREP_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.fasta -out_good ${SEQPREP_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.duplicates_removed -out_bad null -derep 124 -line_width 0
	wait

# Remove FASTA and compress FASTQ files
	rm -f ${SEQPREP_OUTPUT}/${SAMPLE}_all_seqprep.complexity_filtered.fasta
	wait
done
wait
echo "...File setup for MEGAN and/or MIA is complete." >> ${PREFIX}_progress_file_${DATE}.txt

### MIA --- Option 1 -- creating a consensus using iterative mapping - use this with shotgun data
 
#for SAMPLE in $(cat ${PREFIX}-sample-list-${DATE}.txt)
#do
# 	/soe/pheintzman/bin/mia-1.0/src/mia -r ${MIA_REFERENCE_SEQUENCE} -f ${SEQPREP_OUTPUT}/${SAMPLE}_merged.complexity_filtered.fastq -c -C -U -i -F -k 14 -m ${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln >> ${MIA_OUTPUT}/${SAMPLE}_MIA.log.txt 2>&1
 	# wait
 	# gzip ${SEQPREP_OUTPUT}/${SAMPLE}_merged.complexity_filtered.fastq
#wait
#done
 
for SAMPLE in $(cat ${PREFIX}-sample-list-${DATE}.txt)
do 	
	/soe/pheintzman/bin/mia-1.0/src/ma -M ${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln.* -f 3 > ${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln.F.mia_stats.txt
 	wait
 	/soe/pheintzman/bin/mia-1.0/src/ma -M ${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln.* -f 2 > ${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln.F.mia_coverage_per_site.txt
 	wait
 	/soe/pheintzman/bin/mia-1.0/src/ma -M ${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln.* -f 5 > ${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln.F.mia_consensus.fasta
 	wait
 	/soe/pheintzman/bin/mia-1.0/src/ma -M ${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln.* -f 41 > ${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln.F.inputfornext.txt
 	wait
 	perl ${MIA_COVERAGE_FILTER_ANDRE} -c 3 -p 0.67 -I ${SAMPLE}_3x_0.67 <${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln.F.inputfornext.txt >${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln.F.mia_consensus.3x_0.67_filtered.fasta 
  	wait
 	perl ${MIA_COVERAGE_FILTER_ANDRE} -c 10 -p 0.9 -I ${SAMPLE}_10x_0.9 <${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln.F.inputfornext.txt >${MIA_OUTPUT}/${SAMPLE}_merged.complexity_filtered.${MIA_REFERENCE_NAME}.maln.F.mia_consensus.10x_0.9_filtered.fasta 
done
wait
echo "...MIA analyses are complete" >> ${PREFIX}_progress_file_${DATE}.txt 

# cd ${PROCESSING_OUTPUT}

# gzip ${SEQPREP_OUTPUT}/${SAMPLE}_*.fastq
# gzip ${SEQPREP_OUTPUT}/${SAMPLE}_*.fasta
# gzip ${MEGAN_OUTPUT}/${SAMPLE}_all_seqprep.duplicates_removed.BLAST.txt

# echo "...Large files are gzipped."

# Rscript ${CALC_STATS} BWA_analyses/ SeqPrep_output/ MapDamage_output/ MIA_analyses/

# echo "...Summary statistics are estimated!"

# mv ${PREFIX}-sample-list-${DATE}.txt ${PROCESSING_OUTPUT}/Sample_lists_and_progress_files
# mv ${PREFIX}_progress_file_${DATE}.txt ${PROCESSING_OUTPUT}/Sample_lists_and_progress_files

echo "Mitogenome alignment steps are complete."
echo "#######################################################"
