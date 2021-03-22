#!/bin/bash

spack load -r python@3.7.0^gcc@6.3.0
spack load samtools@1.9%gcc@6.3.0
spack load bwa@0.7.15%gcc@6.3.0

# runs all wgs analysis scripts:
# READ IN ARGUMENTS
# usage:
function usage {
        echo "Usage: $(basename $0) [-R] [-G] [-I] [-N]" 2>&1
        echo '   -R   path to reference_genome.fa'
        echo '   -G   path to reference_genome.gbk'
        echo '   -I   path to unmerged input fastq directory'
        echo '   -N   number of processors available'
        exit 1
}

if [[ ${#} -eq 0 ]]; then
   usage
fi

# Define list of arguments expected in the input
optstring="R:G:I:N:"

while getopts ${optstring} arg; do
  case "${arg}" in
    R) REF_FA="${OPTARG}" ;;
    G) REF_GB="${OPTARG}" ;;
    I) IN_DIR="${OPTARG}" ;;
    N) NUM_PROCS="${OPTARG}" ;;

    ?)
      echo "Invalid option: -${OPTARG}."
      echo
      usage
      ;;
  esac
done

## run data processing scripts
source ~/wgs/wgs_scripts_dev/bash_scripts/piecemeal/merge_fastq.sh -I "$IN_DIR"
source ~/wgs/wgs_scripts_dev/bash_scripts/piecemeal/qc_and_trim.sh
source ~/wgs/wgs_scripts_dev/bash_scripts/piecemeal/align_and_clean.sh -R "$REF_FA"

# run SNV callers
source ~/wgs/wgs_scripts_dev/bash_scripts/piecemeal/call_snv.sh -R "$REF_FA"
source ~/wgs/wgs_scripts_dev/bash_scripts/piecemeal/run_snippy.sh -R "$REF_GB" -N "$NUM_PROCS"


# run SNV filters


# run SNV output parsing


# run SV/CNV callers

