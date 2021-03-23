#!/bin/bash

# DEPENDENCIES:
spack load /cjzfz7f  # loads python-2.7.15

function usage {
        echo "Usage: $(basename $0) [-R] -[N]" 2>&1
        echo '   -R   path to reference_genome.fa'
        echo '   -N   number of cores available'
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
mkdir manta
mkdir pilon

# run Manta
echo "%%%%%%%%%%    MANTA     %%%%%%%%%%%"
for f in bam/*dedup.bam
	do
	  BASE=$(basename ${f})
	  mkdir manta/"${BASE/.dedup.bam/}"

    ~biotools/manta-1.6.0.centos6_x86_64/bin/configManta.py \
        --bam "$f" \
        --referenceFasta "$REF" \
        --runDir manta/"${BASE/.dedup.bam/}"

    manta/"${BASE/.dedup.bam/}"/runWorkflow.py -j "$NUM_CORES"

  done


# run Pilon
echo "%%%%%%%%%%    PILON     %%%%%%%%%%%"
for f in bam/*dedup.bam
  do

    BASE=$(basename ${f})
    mkdir pilon/"${BASE/.dedup.bam/}"

    java -jar ~/biotools/pilon-1.24.jar \
      --genome "$REF" --bam "$f" \
      --outdir pilon/"${BASE/.dedup.bam/}" \
      --variant --tracks --nostrays

  done

# run CNGpro
