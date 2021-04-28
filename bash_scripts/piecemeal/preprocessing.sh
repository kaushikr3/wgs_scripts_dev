#!/bin/bash

# DEPENDENCIES:
# fastqc
# trimgalore
# cutadapt
# samtools
spack load fastqc
spack load -r py-cutadapt@1.13
spack load trimgalore

# READ IN ARGUMENTS
# usage:
function usage {
        echo "Usage: $(basename $0) [-I]" 2>&1
        echo '   -I   path to data source directory; holds unmerged_fastq/ dir'
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
mkdir "$SOURCE_DIR"/merged_untrimmed
mkdir "$SOURCE_DIR"/trimmed_fastq
mkdir "$SOURCE_DIR"/reports
mkdir "$SOURCE_DIR"/reports/fastqc_untrimmed_out
mkdir "$SOURCE_DIR"/reports/fastqc_trimmed_out

echo " ----  Processing fastq.gz files from ${SOURCE_DIR}  ----"

# MERGE LANE FASTQS
#for f in "$SOURCE_DIR"/unmerged_fastq/*L001_R1_001.fastq.gz
#do 
#
#      echo "$f"
#      BASE="${f/L001_R1_001.fastq.gz/}"
#      OUTNAME=$(basename ${BASE})
#
#  		FILES_R1=(${BASE}*R1_001.fastq.gz)
#  		FILES_R2=(${BASE}*R2_001.fastq.gz)
#		
#  		echo "Merged ${FILES_R1} & ${FILES_R2} into ${OUTNAME}_R1_001.fastq.gz & ${OUTNAME}_R2_001.fastq.gz"
#		
#  		cat "${FILES_R1[@]}" > ${SOURCE_DIR}/merged_untrimmed/"${OUTNAME}R1_001.fastq.gz" &
#  		cat "${FILES_R2[@]}" > ${SOURCE_DIR}/merged_untrimmed/"${OUTNAME}R2_001.fastq.gz"
#
#done
#
#echo " ---- fastq.gz files merged ----  "
#

# UNTRIMMED FASTQC
#for f in "$SOURCE_DIR"/merged_untrimmed/*_001.fastq.gz
#
#	do
# 			echo Running fastqc on "$f"
#			fastqc -t 4 --extract -o "$SOURCE_DIR"/reports/fastqc_untrimmed_out "$f"
#	done
#

# TRIM AND GENERATE TRIMMED FASTQC
for f in "$SOURCE_DIR"/merged_untrimmed/*R1_001.fastq.gz
	do
			echo "Running TrimGalore on:  ${f}  ${f/R1/R2}"
			trim_galore --paired --output_dir "$SOURCE_DIR"/trimmed_fastq2 \
					"${f}" "${f/R1/R2}"
					#--fastqc --fastqc_args "-t 4 --extract -o "$SOURCE_DIR"/reports/fastqc_trimmed_out" \
	done


spack unload fastqc
spack unload -r py-cutadapt@1.13
spack unload trimgalore



