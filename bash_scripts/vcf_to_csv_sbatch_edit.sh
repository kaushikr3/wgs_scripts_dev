#!/bin/bash

# DEPENDENCIES:
# blast
# conda packages: biopython, pandas, numpy, pyranges

# deal iwth weird python path env bullshit b/c of spack
source ~/.bashrc
conda activate snippy

python_path=$(dirname $(which python))
spack load blast-plus@2.10.0

export PATH="${python_path}:${PATH}"

# READ IN ARGUMENTS
# usage:
function usage {
        echo "Usage: $(basename $0) [-R] [-A] [-L]" 2>&1
        echo '   -R   path to reference_genome.fa'
        echo '   -A   path to annotation_genome.gbk'
        echo '   -L   path to lab references CSV dir'
        exit 1
}

if [[ ${#} -eq 0 ]]; then
   usage
fi

# Define list of arguments expected in the input
optstring="R:A:L:"

while getopts ${optstring} arg; do
  case "${arg}" in
    R) REF="${OPTARG}" ;;
    A) ANNO="${OPTARG}" ;;
    L) LAB="${OPTARG}" ;;

    ?)
      echo "Invalid option: -${OPTARG}."
      echo
      usage
      ;;
  esac
done

# DIRECTORY SETUP
mkdir csv_old
mkdir blast

for f in gatk/*haploid.vcf
	do
 			echo " ----- Running vcf_to_csv on ${f} and ${f/haploid/diploid}----------"

			python ~/wgs/wgs_scripts_dev/python_analysis_scripts/vcf_to_csv.py \
					-vcf "$f" \
					-ref_genome "$REF" -anno_genome "$ANNO" \
					-lab_ref_csv_dir "$LAB" \
					-csv_main_dir csv_old -blast_main_dir_name blast
			
			python ~/wgs/wgs_scripts_dev/python_analysis_scripts/vcf_to_csv.py \
					-vcf "${f/haploid/diploid}" \
					-ref_genome "$REF" -anno_genome "$ANNO" \
					-lab_ref_csv_dir "$LAB" \
					-csv_main_dir csv_old -blast_main_dir_name blast
			
	done


