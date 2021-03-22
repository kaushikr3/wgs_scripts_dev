#!/bin/bash

# DEPENDENCIES:
# bwa
# samtools
# picard

# READ IN ARGUMENTS
# usage:
function usage {
        echo "Usage: $(basename $0) [-R]" 2>&1
        echo '   -R   path to reference_genome.fa'
        exit 1
}

if [[ ${#} -eq 0 ]]; then
   usage
fi

# Define list of arguments expected in the input
optstring="R:"

while getopts ${optstring} arg; do
  case "${arg}" in
    R) REF="${OPTARG}" ;;

    ?)
      echo "Invalid option: -${OPTARG}."
      echo
      usage
      ;;
  esac
done

# DIRECTORY SETUP
mkdir bam

echo " ----  Aligning files against ref: ${REF} ----"


# ALIGNMENT:
for f in fastq/*R1_001_val_1.fq.gz
	do
			# Use MC2155 for WGS of Msmeg
			# Use H37RvCO4 for WGS of H37Rv Mtb to get closest reference genome
			echo "Aligning ${f} and ${f/R1_001_val_1/R2_001_val_2} with BWA MEM" 

			OUTNAME=$(basename ${f})

			bwa mem -M -t 8 "$REF" "$f" "${f/R1_001_val_1/R2_001_val_2}" | \
			    samtools view -Shb - > bam/"${OUTNAME/R1_001_val_1.fq.gz/unsorted.bam}"

	done


# PICARD TOOLS
for f in bam/*unsorted.bam
	do
			echo "Cleaning and sorting ${f} " 
			BASE=$(basename ${f})
			
			# clean unsorted bamfile:
			picard CleanSam \
					INPUT="$f" \
					OUTPUT="${f/unsorted.bam/unsorted.cleaned.bam}"
	
			# sort bam 
			samtools sort -o "${f/unsorted/sorted}" -O bam -T temp -@ 8 "${f/unsorted.bam/unsorted.cleaned.bam}"
		   
			echo "Running MarkDuplicates ${f}" 
			picard MarkDuplicates \
					VALIDATION_STRINGENCY=STRICT \
					INPUT="${f/unsorted/sorted}" \
					OUTPUT="${f/unsorted/strict.dedup}" \
					METRICS_FILE=reports/"${BASE/unsorted.bam/strict.dedup.metrics.txt}" \
					ASSUME_SORTED=true
	
	    	  # removing the add or replace read groups steps here
	#				java -jar ~/biotools/picard/build/libs/picard.jar AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT \
	#						I="$f" O="${f/sorted/RG}" RGID="${f/_sorted.bam/}" RGLB="${a/_sorted.bam/}" RGPL=ILLUMINA \
	#						RGPU="${f/_sorted.bam/}" RGSM="${f/_sorted.bam/}" USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
	
					#java -jar $picard_path MarkDuplicates \
					#java -jar ~/biotools/picard/build/libs/picard.jar MarkDuplicates \
	
			picard MarkDuplicates \
					VALIDATION_STRINGENCY=LENIENT \
					INPUT="${f/unsorted/sorted}" \
					OUTPUT="${f/unsorted/lenient.dedup}" \
					METRICS_FILE=reports/"${BASE/sorted.bam/lenient.dedup.metrics.txt}" \
					ASSUME_SORTED=true
	#				REMOVE_DUPLICATES=TRUE \
	#				USE_JDK_DEFLATER=true \
	#				USE_JDK_INFLATER=true

			echo "Indexing ${f/unsorted/strict.dedup}}"
			samtools index "${f/unsorted/strict.dedup}"

			echo "Indexing ${f/unsorted/lenient.dedup}}"
			samtools index "${f/unsorted/lenient.dedup}"

	done


# EVAL BAMS
for f in bam/*dedup.bam
	do
			echo "BAM  ${f}" >> reports/samtools_stats.log
			samtools flagstat "$f" >> reports/samtools_stats.log
			echo %%%%%%%%%%%%%%% >> reports/samtools_stats.log
			echo  >> reports/samtools_stats.log
	done

