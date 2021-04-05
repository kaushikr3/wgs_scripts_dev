#!/bin/bash

spack load -r python@3.7.0^gcc@6.3.0
spack load samtools@1.9%gcc@6.3.0
spack load bwa@0.7.15%gcc@6.3.0

# runs all wgs analysis scripts:
# READ IN ARGUMENTS
# usage:
function usage {
        echo "Usage: $(basename $0) [-R] [-I] [-N]" 2>&1
        echo '   -R   path to reference_genome.fa'
        echo '   -I   path to unmerged input fastq directory'
        echo '   -N   number of processors available'
        exit 1
}

if [[ ${#} -eq 0 ]]; then
   usage
fi

# Define list of arguments expected in the input
optstring="R:I:N:"

while getopts ${optstring} arg; do
  case "${arg}" in
    R) REF="${OPTARG}" ;;
    I) IN_DIR="${OPTARG}" ;;
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
#source ~/wgs/wgs_scripts_dev/bash_scripts/piecemeal/merge_fastq.sh -I "$IN_DIR"
#source ~/wgs/wgs_scripts_dev/bash_scripts/piecemeal/qc_and_trim.sh
source ~/wgs/wgs_scripts_dev/bash_scripts/piecemeal/align_and_clean.sh -R "$REF" -N "$NUM_CORES"
#
## run SNV callers
source ~/wgs/wgs_scripts_dev/bash_scripts/piecemeal/call_snv.sh -R "$REF"
source ~/wgs/wgs_scripts_dev/bash_scripts/piecemeal/run_snippy.sh -R "$REF" -N "$NUM_CORES"

# run SNV filters
#source ~/wgs/wgs_scripts_dev/bash_scripts/piecemeal/call_structural_variants.sh -R "$REF"

# run SNV output parsing
# run SV/CNV callers


