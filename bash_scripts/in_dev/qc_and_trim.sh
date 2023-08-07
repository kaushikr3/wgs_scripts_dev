#!/bin/bash

shopt -s nullglob

# DEPENDENCIES:
# fastqc
# trimgalore
# cutadapt
# samtools
spack load fastqc
spack load -r py-cutadapt@1.13
spack load trimgalore


# DIRECTORY SETUP
#mkdir fastq
#mkdir fastq/merged_untrimmed

mkdir reports
mkdir reports/fastqc_untrimmed_out
mkdir reports/fastqc_trimmed_out

echo " ----  Trimming and Running QC  ----"

# UNTRIMMED FASTQC
#for f in fastq/merged_untrimmed/*_001.fastq.gz
#
#	do
# 			echo Running fastqc on "$f"
#			fastqc -t 2 --extract -o reports/fastqc_untrimmed_out "$f"
#	done
#

# TRIM AND GENERATE TRIMMED FASTQC
#for f in fastq/merged_untrimmed/*R1_001.fastq.gz
for f in raw_fastq/*R1.fastq.gz
	do
			echo "Running TrimGalore on:  ${f}  ${f/R1/R2}"
			trim_galore --paired --output_dir trimmed_fastq \
					"${f}" "${f/R1/R2}"
					#--fastqc --fastqc_args "-t 2 --extract -o reports/fastqc_trimmed_out" \
	done


spack unload fastqc
spack unload -r py-cutadapt@1.13
spack unload trimgalore


