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
mkdir freebayes

echo " ----  Calling SNVs with reference: ${REF} ----"

# FREEBAYES

for f in bam/*dedup.bam
do
	   	echo "Running freebayes on ${f}" 
	   	BASE=$(basename "${f}")

	   	~/biotools/freebayes-1.3.4-linux-static-AMD64 -f "$REF" --ploidy 1 "$f" | bcftools filter \
				-O v -e 'DP<=5' -e 'SAF == 0 || SAR == 0' \
				-o freebayes/"${BASE/dedup.bam/freebayes.filt.vcf}" -

		
done

spack unload bcftools@1.9%gcc@6.3.0


