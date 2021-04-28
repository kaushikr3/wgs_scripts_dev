#!/bin/bash

# DEPENDENCIES:
# none

# READ IN ARGUMENTS
# usage:
function usage {
        echo "Usage: $(basename $0) [-I]" 2>&1
        echo '   -I   path to unmerged_fastq directory: unmerged_fastq/*_unmerged.fastq'
        exit 1
}

if [[ ${#} -eq 0 ]]; then
   usage
fi

# Define list of arguments expected in the input
optstring="I:"

while getopts ${optstring} arg; do
  case "${arg}" in
    I) SOURCE_DIR="${OPTARG}" ;;

    ?)
      echo "Invalid option: -${OPTARG}."
      echo
      usage
      ;;
  esac
done

# DIRECTORY SETUP
#mkdir fastq
#mkdir fastq/merged_untrimmed

echo " ----  Processing fastq.gz files from ${SOURCE_DIR}  ----"

# MERGE LANE FASTQS
for f in ${SOURCE_DIR}/unmerged_fastq/*L001_R1_001.fastq.gz
do 

      echo "$f"
      BASE="${f/L001_R1_001.fastq.gz/}"
      OUTNAME=$(basename ${BASE})

  		FILES_R1=(${BASE}*R1_001.fastq.gz)
  		FILES_R2=(${BASE}*R2_001.fastq.gz)
		
  		echo "Merged ${FILES_R1} & ${FILES_R2} into ${OUTNAME}_R1_001.fastq.gz & ${OUTNAME}_R2_001.fastq.gz"
		
  		cat "${FILES_R1[@]}" > ${SOURCE_DIR}/merged_untrimmed/"${OUTNAME}R1_001.fastq.gz" &
  		cat "${FILES_R2[@]}" > ${SOURCE_DIR}/merged_untrimmed/"${OUTNAME}R2_001.fastq.gz"

done

echo " ---- fastq.gz files merged ----  "
