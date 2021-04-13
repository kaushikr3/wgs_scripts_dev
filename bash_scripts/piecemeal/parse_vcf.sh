#!/bin/bash

# DEPENDENCIES:
# vcflib
# gatk
spack load miniconda3@4.6.14
conda activate snippy


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

echo "Left aligning and parsing VCFs"

for f in vcf/*.vcf
do

	   	vcfleftalign -r "$REF" "$f" > "${f/vcf/leftalign.vcf}"

      ~/biotools/gatk-4.2.0.0/gatk VariantsToTable -V "${f/vcf/leftalign.vcf}" \
        -F CHROM -F POS -F REF -F ALT -F QUAL -F AC -F AF -F QD \
        -GF AD -GF DP -GF GQ -GF GT -GF PL \
        -O "${f/vcf/tsv}"

done

conda deactivate
spack unload miniconda3@4.6.14
