#!/bin/bash

spack load -r python@3.7.0^gcc@6.3.0
spack load samtools@1.9%gcc@6.3.0
spack load bwa@0.7.15%gcc@6.3.0

# runs all wgs analysis scripts:
# READ IN ARGUMENTS
# usage:
function usage {
        echo "Usage: $(basename $0) [-R] [-I] [-F] [-N]" 2>&1
        echo '   -R   path to reference_genome.fa'
        echo '   -I   path to unmerged input fastq directory'
		echo '   -F   path to trimmed fastq directory (optional)'
        echo '   -N   number of processors available'
        exit 1
}

if [[ ${#} -eq 0 ]]; then
   usage
fi

# Define default FASTQ_DIR:
FASTQ_DIR=fastq

# Define list of arguments expected in the input
optstring="R:I:F:N:"

while getopts ${optstring} arg; do
  case "${arg}" in
    R) REF_FA="${OPTARG}" ;;
    I) IN_DIR="${OPTARG}" ;;
    F) FASTQ_DIR="${OPTARG}" ;;
    N) NUM_CORES="${OPTARG}" ;;

    ?)
      echo "Invalid option: -${OPTARG}."
      echo
      usage
      ;;
  esac
done

echo "$NUM_CORES"

## run data processing scripts
#~/wgs/wgs_scripts_dev/bash_scripts/piecemeal/preprocessing.sh -I "$IN_DIR"
#~/wgs/wgs_scripts_dev/bash_scripts/piecemeal/merge_fastq.sh -I "$IN_DIR"
#~/wgs/wgs_scripts_dev/bash_scripts/piecemeal/qc_and_trim.sh

~/wgs/wgs_scripts_dev/bash_scripts/piecemeal/align_and_clean.sh -R "$REF_FA" -I "$FASTQ_DIR" -N "$NUM_CORES"

## run snippy 
#sbatch --export=R="$REF_FA",G="$REF_GBK" ~/wgs/wgs_scripts_dev/bash_scripts/piecemeal/sbatch_snippy.sh

## run variant callers
#~/wgs/wgs_scripts_dev/bash_scripts/piecemeal/call_snv.sh -R "$REF_FA"
#~/wgs/wgs_scripts_dev/bash_scripts/piecemeal/call_structural_variants.sh -R "$REF_FA"

## run vcf parsing
#~/wgs/wgs_scripts_dev/bash_scripts/piecemeal/left_align_parse_vcf.sh -R "$REF_FA"


