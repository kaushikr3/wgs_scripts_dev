#!/bin/bash

# DEPENDENCIES:
# bwa
# samtools
# picard
# gatk
# freebayes


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
mkdir gatk
mkdir freebayes

echo " ----  Calling SNVs with reference: ${REF} ----"

# GATK ANALYSIS
for f in bam/*dedup.bam
	do
			echo "Running GATK on ${f}" 
			BASE=$(basename ${f})
			
			~/biotools/gatk-4.2.0.0/gatk HaplotypeCaller \
					-R "$REF" -I "$f" \
				   	--sample-ploidy 1 \
				   	-O gatk/"${BASE/dedup.bam/gatk.haploid.vcf}"
		   
			~/biotools/gatk-4.2.0.0/gatk HaplotypeCaller \
				   	-R "$REF" -I "$f" \
				   	-O gatk/"${BASE/dedup.bam/gatk.diploid.vcf}"

	done


# FREEBAYES

for f in bam/*.dedup.bam
do
	   	echo "Running freebayes on ${f}" 
	   	BASE=$(basename "${f}")

	   	~/biotools/freebayes-1.3.4-linux-static-AMD64 -f "$REF" --ploidy 1 "$f" > freebayes/"${BASE/dedup.bam/freebayes.vcf}"
done
