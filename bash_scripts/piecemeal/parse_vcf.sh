#!/bin/bash

# DEPENDENCIES:
# vcflib
# gatk

source ~/.bashrc

#spack load miniconda3@4.6.14
#conda activate snippy


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

#
#echo "Parsing gatk VCFs" 
#
#for f in gatk/*.filt.vcf
#do 
#		BASE=$(basename ${f})	
#		
#		~/biotools/gatk-4.2.0.0/gatk VariantsToTable -V "$f" \
#				-F CHROM -F POS -F REF -F ALT -F QUAL -F AC -F AF -F QD \
#				-GF AD -GF DP -GF GQ -GF GT -GF PL \
#				-O vcf_parse/"${BASE/.vcf/.tsv}"
#
#done
#
#echo "Parsing freebayes VCFs"
#
#for f in freebayes/*.filt.vcf
#do 
#		BASE=$(basename ${f})	
#		
#		~/biotools/gatk-4.2.0.0/gatk VariantsToTable -V "$f" \
#				-F CHROM -F POS -F REF -F ALT -F QUAL -F AC -F AF -F QD \
#				-GF AD -GF DP -GF GQ -GF GT -GF PL \
#				-O vcf_parse/"${BASE/.vcf/.tsv}"
#
#done
#
#
#echo "Parsing pilon VCFs"
#for f in pilon/*_pilon.vcf
#do 
#		BASE=$(basename ${f})	
#		
#		~/biotools/gatk-4.2.0.0/gatk VariantsToTable -V "$f" \
#				-F CHROM -F POS -F REF -F ALT -F QUAL \
#				-F DP -F TD -F BQ -F QD -F IC -F DC -F AC -F AF -F SVTYPE -F SVLEN -F END \
#				-GF GT -GF AD -GF DP \
#				-O vcf_parse/"${BASE/.vcf/.tsv}"
#
#done
#
echo "Parsing snippy VCFs"
for f in snippy/*.filt.leftalign.vcf
do 
		BASE=$(basename ${f})	
		
		~/biotools/gatk-4.2.0.0/gatk VariantsToTable -V "$f" \
				-F CHROM -F POS -F REF -F ALT -F QUAL -F AC -F AF -F QD \
				-GF AD -GF DP -GF GQ -GF GT -GF PL \
				-O leftalign_vcf_parse/"${BASE/.vcf/.tsv}"

done


