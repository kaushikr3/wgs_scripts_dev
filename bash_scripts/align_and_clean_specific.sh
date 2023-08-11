#!/bin/bash

# DEPENDENCIES:
# bwa
# samtools
# picard
spack load picard

# READ IN ARGUMENTS
# usage:
function usage {
        echo "Usage: $(basename $0) [-R] [-I] [-N] [-F] [-1] [-2]" 2>&1
        echo '   -R   path to reference_genome.fa'
        echo '   -I   relative path to trimmed fastq dir'
        echo '   -N   number of cores avaliable'
        echo '   -F   file list to process with R1_001_val_1.fq.gz filenames'
        echo '   -1  string at end of R1 trimmed fastq file to  replace'
        echo '   -2  string at end of R2 trimmed fastq file to  replace'
        exit 1
}

if [[ ${#} -eq 0 ]]; then
   usage
fi

# Default directory, if none input:
IN=fastq
R1=R1_001_val_1.fq.gz
R2=R2_001_val_2.fq.gz

# Define list of arguments expected in the input
optstring="R:I:N:F:1:2:"

while getopts ${optstring} arg; do
  case "${arg}" in
    R) REF="${OPTARG}" ;;
    I) IN="${OPTARG}" ;;
    N) NUM_CORES="${OPTARG}" ;;
	F) sample_list="${OPTARG}" ;;
	1) R1="${OPTARG}" ;;
	2) R2="${OPTARG}" ;;

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

			echo "Aligning ${f} and ${f/"$R1"/$R2} with BWA MEM" 

			bwa mem -M -t "$NUM_CORES" "$REF" "$IN"/"$f" "$IN"/"${f/"$R1"/$R2}" | \
			    samtools view -Shb - > bam/"${f/"$R1"/unsorted.bam}"

			#OUTNAME=$(basename ${f})
			#echo "Aligning ${f} and ${f/R1_001_val_1/2_val_2} with BWA MEM" 
			#bwa mem -M -t "$NUM_CORES" "$REF" "$IN"/"$f" "$IN"/"${f/R1_001_val_1/2_val_2}" | \
			#    samtools view -Shb - > bam/"${f/R1_001_val_1.fq.gz/unsorted.bam}"

	done < "$sample_list"


# PICARD TOOLS
#for f in bam/*unsorted.bam

while read f
	do
			echo "Cleaning and sorting ${f/"$R1"/} " 
			#BASE=$(basename ${f})
			
			# clean unsorted bamfile:
			picard CleanSam \
					INPUT=bam/"${f/"$R1"/unsorted.bam}" \
					OUTPUT=bam/"${f/"$R1"/unsorted.cleaned.bam}"
			
			picard AddOrReplaceReadGroups \
					I=bam/"${f/"$R1"/unsorted.cleaned.bam}" \
					O=bam/"${f/"$R1"/RG.bam}" \
          			RGID="${f/"$R1"/unsorted.bam/}" \
          			RGPL=ILLUMINA \
          			RGLB=lib \
          			RGPU=unit \
         			RGSM=sample

#          			RGSM="${BASE/unsorted.bam/}"
#         			RGLB="${a/_sorted.bam/}" \
#         			RGPU="${f/_sorted.bam/}" \
#         			USE_JDK_DEFLATER=true USE_JDK_INFLATER=true

			# sort bam 
			samtools sort -o bam/"${f/"$R1"/sorted.bam}" -O bam -T "${f/"$R1"/}" -@ 8 bam/"${f/"$R1"/RG.bam}"
			#samtools sort -o bam/"${f/R1_001_val_1.fq.gz/sorted.bam}" -O bam -T temp -@ 8 bam/"${f/R1_001_val_1.fq.gz/RG.bam}"
		   
			echo "Running MarkDuplicates ${f}" 
			picard MarkDuplicates \
					INPUT=bam/"${f/"$R1"/sorted.bam}" \
					OUTPUT=bam/"${f/"$R1"/dedup.bam}" \
					METRICS_FILE=reports/"${f/"$R1"/dedup.metrics.txt}" \
					ASSUME_SORTED=true
	
			# index bam 
			samtools index bam/"${f/"$R1"/dedup.bam}" 

	done < "$sample_list"


# EVAL BAMS
#for f in bam/*dedup.bam

while read f
	do

			echo "Indexing ##################"
			echo "BAM  bam/${f/"$R1"/dedup.bam}" >> reports/samtools_stats.log
			
			#samtools index bam/"${f/R1_001_val_1.fq.gz/dedup.bam}" 
			samtools flagstat bam/"${f/"$R1"/dedup.bam}" >> reports/samtools_stats.log

			echo %%%%%%%%%%%%%%% >> reports/samtools_stats.log
			echo  >> reports/samtools_stats.log

	done < "$sample_list"


spack unload picard

