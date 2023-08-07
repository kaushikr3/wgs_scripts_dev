#!/bin/bash

# DEPENDENCIES:
# bwa
# samtools
# picard
# gatk
# freebayes
#spack load bcftools@1.9%gcc@6.3.0

# READ IN ARGUMENTS
# usage:
function usage {
        echo "Usage: $(basename $0) [-R] [-F]" 2>&1
        echo '   -R   path to reference_genome.fa'
        echo '   -F   path to file holding sample names of *dedup.bam type'
        exit 1
}

if [[ ${#} -eq 0 ]]; then
   usage
fi

# Define list of arguments expected in the input
optstring="R:F:"

while getopts ${optstring} arg; do
  case "${arg}" in
    R) REF="${OPTARG}" ;;
    F) file_list="${OPTARG}" ;;

    ?)
      echo "Invalid option: -${OPTARG}."
      echo
      usage
      ;;
  esac
done

# DIRECTORY SETUP
mkdir gatk

echo " ----  Calling SNVs with reference: ${REF} ----"

# GATK ANALYSIS
while read f
	do
			echo "Running GATK on ${f}" 
			
			~/biotools/gatk-4.2.0.0/gatk HaplotypeCaller --sample-ploidy 1 \
					-R "$REF" -I bam/"$f" \
				   	-O gatk/"${f/dedup.bam/gatk.haploid.vcf}"
		   
			~/biotools/gatk-4.2.0.0/gatk HaplotypeCaller \
				   	-R "$REF" -I bam/"$f" \
				   	-O gatk/"${f/dedup.bam/gatk.diploid.vcf}"

	done < $file_list



