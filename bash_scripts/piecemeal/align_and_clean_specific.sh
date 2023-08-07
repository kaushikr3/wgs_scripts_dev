#!/bin/bash

# DEPENDENCIES:
# bwa
# samtools
# picard
spack load picard

# READ IN ARGUMENTS
# usage:
function usage {
        echo "Usage: $(basename $0) [-R] [-I] [-N] [-F]" 2>&1
        echo '   -R   path to reference_genome.fa'
        echo '   -I   relative path to trimmed fastq dir'
        echo '   -N   number of cores avaliable'
        echo '   -F   file list to process with 1_val_1.fq.gz filenames'
        exit 1
}

if [[ ${#} -eq 0 ]]; then
   usage
fi

# Default directory, if none input:
IN=fastq

# Define list of arguments expected in the input
optstring="R:I:N:F:"

while getopts ${optstring} arg; do
  case "${arg}" in
    R) REF="${OPTARG}" ;;
    I) IN="${OPTARG}" ;;
    N) NUM_CORES="${OPTARG}" ;;
	F) file_list="${OPTARG}" ;;

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

while read f
	do

			#echo "Aligning ${f} and ${f/R1_001_val_1/R2_001_val_2} with BWA MEM" 
			echo "Aligning ${f} and ${f/1_val_1/2_val_2} with BWA MEM" 
			#OUTNAME=$(basename ${f})

			#bwa mem -M -t "$NUM_CORES" "$REF" "$f" "${f/R1_001_val_1/R2_001_val_2}" | \
			#    samtools view -Shb - > bam/"${OUTNAME/R1_001_val_1.fq.gz/unsorted.bam}"

			bwa mem -M -t "$NUM_CORES" "$REF" "$IN"/"$f" "$IN"/"${f/1_val_1/2_val_2}" | \
			    samtools view -Shb - > bam/"${f/1_val_1.fq.gz/unsorted.bam}"

	done < "$file_list"


# PICARD TOOLS
#for f in bam/*unsorted.bam

while read f
	do
			echo "Cleaning and sorting ${f} " 
			#BASE=$(basename ${f})
			
			# clean unsorted bamfile:
			picard CleanSam \
					INPUT=bam/"${f/1_val_1.fq.gz/unsorted.bam}" \
					OUTPUT=bam/"${f/1_val_1.fq.gz/unsorted.cleaned.bam}"
			
			picard AddOrReplaceReadGroups \
					I=bam/"${f/1_val_1.fq.gz/unsorted.cleaned.bam}" \
					O=bam/"${f/1_val_1.fq.gz/RG.bam}" \
          			RGID="${f/1_val_1.fq.gz/unsorted.bam/}" \
          			RGPL=ILLUMINA \
          			RGLB=lib \
          			RGPU=unit \
         			RGSM=sample

#          			RGSM="${BASE/unsorted.bam/}"
#         			RGLB="${a/_sorted.bam/}" \
#         			RGPU="${f/_sorted.bam/}" \
#         			USE_JDK_DEFLATER=true USE_JDK_INFLATER=true

			# sort bam 
			samtools sort -o bam/"${f/1_val_1.fq.gz/sorted.bam}" -O bam -T "${f/1_val_1.fq.gz/}" -@ 8 bam/"${f/1_val_1.fq.gz/RG.bam}"
			#samtools sort -o bam/"${f/1_val_1.fq.gz/sorted.bam}" -O bam -T temp -@ 8 bam/"${f/1_val_1.fq.gz/RG.bam}"
		   
			echo "Running MarkDuplicates ${f}" 
			picard MarkDuplicates \
					INPUT=bam/"${f/1_val_1.fq.gz/sorted.bam}" \
					OUTPUT=bam/"${f/1_val_1.fq.gz/dedup.bam}" \
					METRICS_FILE=reports/"${f/1_val_1.fq.gz/dedup.metrics.txt}" \
					ASSUME_SORTED=true
	
			# index bam 
			samtools index bam/"${f/1_val_1.fq.gz/dedup.bam}" 

	done < "$file_list"


# EVAL BAMS
#for f in bam/*dedup.bam

while read f
	do

			echo "Indexing ##################"
			echo "BAM  bam/${f/1_val_1.fq.gz/dedup.bam}" >> reports/samtools_stats.log
			
			#samtools index bam/"${f/1_val_1.fq.gz/dedup.bam}" 
			samtools flagstat bam/"${f/1_val_1.fq.gz/dedup.bam}" >> reports/samtools_stats.log

			echo %%%%%%%%%%%%%%% >> reports/samtools_stats.log
			echo  >> reports/samtools_stats.log

	done < "$file_list"


spack unload picard
