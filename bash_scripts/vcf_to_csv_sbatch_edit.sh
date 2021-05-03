#!/bin/bash

# DEPENDENCIES:
# blast
# conda packages: biopython, pandas, numpy, pyranges

source ~/.bashrc
conda activate snippy

python_path=$(dirname $(which python))
echo "PYTHON:"
echo "$python_path"
python --version

spack load blast-plus@2.10.0
echo "BLAST CALLED"

echo "PATH:"
export PATH="${python_path}:${PATH}"
echo "$PATH"

echo "PYTHON:"
which python
python --version



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

for f in gatk_copies/*haploid.vcf
	do
 			echo " ----- Running vcf_to_csv on ${f} and ${f/haploid/diploid}----------"

			python ~/wgs/wgs_scripts_dev/python_analysis_scripts/vcf_to_csv.py \
					-vcf "$f" \
					-ref_genome "$REF" -anno_genome "$ANNO" \
					-lab_ref_csv_dir "$LAB" \
					-csv_main_dir ~/csv_temp_dl -blast_main_dir_name ~/blast_temp_dl
					#-csv_main_dir ~/csv_temp_dl -blast_main_dir_name ~/blast_temp_dl
			
			python ~/wgs/wgs_scripts_dev/python_analysis_scripts/vcf_to_csv.py \
					-vcf "${f/haploid/diploid}" \
					-ref_genome "$REF" -anno_genome "$ANNO" \
					-lab_ref_csv_dir "$LAB" \
					-csv_main_dir ~/csv_temp_dl -blast_main_dir_name ~/blast_temp_dl
					#-csv_main_dir ~/csv_temp -blast_main_dir_name ~/blast_temp
			
	done


