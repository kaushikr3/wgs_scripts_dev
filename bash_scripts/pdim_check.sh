#!/bin/bash

# DEPENDENCIES:
# conda packages: biopython, pandas, numpy, pyranges

source ~/.bashrc
conda activate snippy

# READ IN ARGUMENTS
# usage:
function usage {
        echo "Usage: $(basename $0) [-R] [-S] [-C]" 2>&1
        echo '   -R   path to reference_genome.gbk'
        echo '   -S   Background strain: H37RvCO, HN878, Erdman, BCG'
        echo '   -C   Csv directory path'
        exit 1
}

if [[ ${#} -eq 0 ]]; then
   usage
fi

# set default csv_dir value:
CSV_DIR='snp_xlsx'
#CSV_DIR='csv'

# Define list of arguments expected in the input
optstring="R:S:C:"

while getopts ${optstring} arg; do
  case "${arg}" in
    R) REF="${OPTARG}" ;;
    S) STRAIN="${OPTARG}" ;;
    C) CSV_DIR="${OPTARG}" ;;

    ?)
      echo "Invalid option: -${OPTARG}."
      echo
      usage
      ;;
  esac
done

# DIRECTORY SETUP
mkdir cov

echo "Generating Coverage Files"

#for f in bam/*dedup.bam
#	do
#			BASE=$(basename ${f})
#			echo "$BASE"
#			bedtools genomecov -ibam "$f" -d > cov/${BASE/_dedup.bam/.cov}
#	done
#
echo "Running pdim_check"

python ~/wgs/wgs_scripts_dev/python_analysis_scripts/pdim_check_xlsx.py \
		-strain "$STRAIN" -stringency Stringent -csv_dir "$CSV_DIR" -genome "$REF"

#python ~/wgs/wgs_scripts_dev/python_analysis_scripts/pdim_check.py \
#		-strain "$STRAIN" -stringency stringent -csv_dir "$CSV_DIR" -genome "$REF"
#


