#!/bin/bash

# DEPENDENCIES:
# ???

# READ IN ARGUMENTS
# usage:
function usage {
        echo "Usage: $(basename $0) [-R]" 2>&1
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
    N) NUM_PROCS="${OPTARG}" ;;

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

# GATK ANALYSIS
for f in fastq/*R1_001_val_1.fq.gz
	do

			OUTNAME=$(basename ${f})

      snippy --cpus "$NUM_PROCS" --outdir snippy --ref "$REF" --prefix $OUTNAME --cleanup\
          --R1 "$f" --R2 "${f/R1_001_val_1/R2_001_val_2}"

	done
