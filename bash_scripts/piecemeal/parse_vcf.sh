#!/bin/bash

# DEPENDENCIES:
# vcflib
# gatk

source ~/.bashrc

#spack load miniconda3@4.6.14
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


mkdir vcf_parse
#mkdir vcf_parse/gatk
#mkdir vcf_parse/freebayes

echo "Parsing gatk VCFs" 

for f in gatk/*gatk.haploid.vcf
do 
		BASE=$(basename ${f})	

		#vcfleftalign -r "$REF" "$f" > "${f/.vcf/.leftalign.vcf}"
		
		#~/biotools/gatk-4.2.0.0/gatk VariantsToTable -V "${f/.vcf/.leftalign.vcf}" \
		~/biotools/gatk-4.2.0.0/gatk VariantsToTable -V "$f" \
				-F CHROM -F POS -F REF -F ALT -F QUAL \
				-GF GT -GF GQ -GF PL -GF AD -GF DP \
				-O vcf_parse/"${BASE/.vcf/.tsv}"

done

for f in gatk/*gatk.diploid.vcf
do 
		BASE=$(basename ${f})	

		#vcfleftalign -r "$REF" "$f" > "${f/.vcf/.leftalign.vcf}"
		
		#~/biotools/gatk-4.2.0.0/gatk VariantsToTable -V "${f/.vcf/.leftalign.vcf}" \
		~/biotools/gatk-4.2.0.0/gatk VariantsToTable -V "$f" \
				-F CHROM -F POS -F REF -F ALT -F QUAL \
				-GF GT -GF GQ -GF PL -GF AD -GF DP \
				-O vcf_parse/"${BASE/.vcf/.tsv}"

done

#echo "Parsing freebayes VCFs"
#
#for f in freebayes/*freebayes.vcf
#do 
#		BASE=$(basename ${f})	
#
#		vcfleftalign -r "$REF" "$f" > "${f/.vcf/.leftalign.vcf}"
#
#		bcftools +tag2tag "${f/.vcf/.leftalign.vcf}" -- -r --gl-to-pl > "${f/vcf/PL.vcf}"
#		~/biotools/gatk-4.2.0.0/gatk VariantsToTable \
#				-V "${f/vcf/PL.vcf}" \
#				-F CHROM -F POS -F REF -F ALT -F QUAL \
#				-GF GT -GF GQ -GF PL -GF AD -GF DP \
#				-O vcf_parse/freebayes/"${BASE/.vcf/.tsv}"
#
#done
#
#
#mkdir vcf_parse/pilon
#
#echo "Parsing pilon VCFs"
#for f in pilon/*_pilon.vcf
#do 
#		BASE=$(basename ${f})	
#
#		#vcfleftalign -r "$REF" "$f" > "${f/.vcf/.leftalign.vcf}"
#		
#		~/biotools/gatk-4.2.0.0/gatk VariantsToTable -V "$f" \
#				-F CHROM -F POS -F REF -F ALT -F QUAL -F FILTER \
#				-F PC -F BQ -F QD -F IC -F DC -F XC -F AC -F AF -F SVTYPE -F SVLEN -F END -F IMPRECISE \
#				-GF GT -GF AD -GF DP \
#				-O vcf_parse/pilon/"${BASE/.vcf/.tsv}"
#
#				#-V "${f/.vcf/.leftalign.vcf}" \
#
#done


#echo "Parsing snippy VCFs"
#for f in snippy/*.filt.leftalign.vcf
#do 
#		BASE=$(basename ${f})	
#		bcftools +tag2tag "$f" -- -r --gl-to-pl > "${f/vcf/PL.vcf}"
#		~/biotools/gatk-4.2.0.0/gatk VariantsToTable \
#				-V "${f/vcf/PL.vcf}" \
#				-F CHROM -F POS -F REF -F ALT -F QUAL \
#				-GF GT -GF PL -GF DP -GF RO -GF AO \
#				-O vcf_parse/"${BASE/.vcf/.tsv}"
#
#done
