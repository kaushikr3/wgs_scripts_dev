#!/bin/bash

# DEPENDENCIES:
spack load miniconda3@4.6.14
conda activate snippy

# READ IN ARGUMENTS
# usage:
function usage {
        echo "Usage: $(basename $0) [-R] [-N]" 2>&1
        echo '   -R   path to reference_genome.gbk'
        echo '   -N   number of processors available'
        exit 1
}

if [[ ${#} -eq 0 ]]; then
   usage
fi

# Define list of arguments expected in the input
optstring="R:N:"

while getopts ${optstring} arg; do
  case "${arg}" in
    R) REF="${OPTARG}" ;;
    N) NUM_CORES="${OPTARG}" ;;

    ?)
      echo "Invalid option: -${OPTARG}."
      echo
      usage
      ;;
  esac
done

# DIRECTORY SETUP
mkdir snippy

echo " ----  Running SNIPPY with reference: ${REF} ----"

# SNIPPY ANALYSIS
for f in fastq/*R1_001_val_1.fq.gz
	do
		   	BASE=$(basename ${f}) 
			
			snippy --cpus "$NUM_CORES" --cleanup --force \
					--outdir snippy -prefix "${BASE/_R1_001_val_1.fq.gz/}" \
					--reference "$REF" \
					--R1 "$f" --R2 "${f/R1_001_val_1/R2_001_val_2}"

	done


echo "Left aligning and parsing snippy VCFs"
for f in snippy/*.filt.vcf
do
		BASE=$(basename ${f})
		vcfleftalign -r "$REF" "$f" > "${f/.vcf/.leftalign.vcf}"
		
		~/biotools/gatk-4.2.0.0/gatk VariantsToTable -V "${f/.vcf/.leftalign.vcf}" \
				-F CHROM -F POS -F REF -F ALT -F QUAL -F AC -F AF -F QD \
				-GF AD -GF DP -GF GQ -GF GT -GF PL \
				-O vcf_parse/"${BASE/.vcf/_snippy.tsv}"

done



conda deactivate
spack unload miniconda3@4.6.14


