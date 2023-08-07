#!/bin/bash

# DEPENDENCIES:
# bwa
# samtools
# picard
spack load picard

# READ IN ARGUMENTS
# usage:
function usage {
        echo "Usage: $(basename $0) [-R] [-I] [-N]" 2>&1
        echo '   -R   path to reference_genome.fa'
        echo '   -I   relative path to trimmed fastq dir'
        echo '   -N   number of cores avaliable'
        exit 1
}

if [[ ${#} -eq 0 ]]; then
   usage
fi

# Default directory, if none input:
IN=fastq

# Define list of arguments expected in the input
optstring="R:I:N:"

while getopts ${optstring} arg; do
  case "${arg}" in
    R) REF="${OPTARG}" ;;
    I) IN="${OPTARG}" ;;
    N) NUM_CORES="${OPTARG}" ;;

    ?)
      echo "Invalid option: -${OPTARG}."
      echo
      usage
      ;;
  esac
done

# DIRECTORY SETUP
mkdir bam
mkdir reports

echo " ----  Aligning files against ref: ${REF} ----"


# ALIGNMENT:
#for f in "$IN"/*R1_001_val_1.fq.gz
#for f in "$IN"/*1_val_1.fq.gz
#	do
#
#			#echo "Aligning ${f} and ${f/R1_001_val_1/R2_001_val_2} with BWA MEM" 
#			echo "Aligning ${f} and ${f/1_val_1/2_val_2} with BWA MEM" 
#			OUTNAME=$(basename ${f})
#
#			#bwa mem -M -t "$NUM_CORES" "$REF" "$f" "${f/R1_001_val_1/R2_001_val_2}" | \
#			#    samtools view -Shb - > bam/"${OUTNAME/R1_001_val_1.fq.gz/unsorted.bam}"
#
#			bwa mem -M -t "$NUM_CORES" "$REF" "$f" "${f/1_val_1/2_val_2}" | \
#			    samtools view -Shb - > bam/"${OUTNAME/1_val_1.fq.gz/unsorted.bam}"
#
#	done


# PICARD TOOLS
for f in bam2/*unsorted.bam
	do
			echo "Cleaning and sorting ${f} " 
			BASE=$(basename ${f})
			
			# clean unsorted bamfile:
			picard CleanSam \
					INPUT="$f" \
					OUTPUT="${f/unsorted/unsorted.cleaned}" 
			
			picard AddOrReplaceReadGroups \
          			I="${f/unsorted/unsorted.cleaned}" O="${f/unsorted/RG}" \
          			RGID="${BASE/unsorted.bam/}" \
          			RGPL=ILLUMINA \
          			RGLB=lib \
          			RGPU=unit \
         			RGSM=sample

#          			RGSM="${BASE/unsorted.bam/}"
#         			RGLB="${a/_sorted.bam/}" \
#         			RGPU="${f/_sorted.bam/}" \
#         			USE_JDK_DEFLATER=true USE_JDK_INFLATER=true

			# sort bam 
			samtools sort -o "${f/unsorted/sorted}" -O bam -T temp -@ 8 "${f/unsorted/RG}"
		   
			echo "Running MarkDuplicates ${f}" 
			picard MarkDuplicates \
					INPUT="${f/unsorted/sorted}" \
					OUTPUT="${f/unsorted/dedup}" \
					METRICS_FILE=reports/"${BASE/unsorted.bam/dedup.metrics.txt}" \
					ASSUME_SORTED=true
	

	done


# INDEX AND EVAL BAMS
for f in bam2/*dedup.bam
	do

			echo "Indexing ##################"
			echo "BAM  ${f}" >> reports/samtools_stats.log
			
			samtools index "${f/unsorted/dedup}" &
			samtools flagstat "$f" >> reports/samtools_stats.log

			echo %%%%%%%%%%%%%%% >> reports/samtools_stats.log
			echo  >> reports/samtools_stats.log

	done


spack unload picard
