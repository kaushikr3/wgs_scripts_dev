#!/bin/bash

#mkdir fastqc_untrimmed_out
#mkdir merged_fastq
mkdir trimmed_fastq

# merge fastq files
for a in $(ls individual_fastq/*L001_R1_001.fastq.gz)
	do
		BASENAME=$(basename ${a})
		cat ${a} ${a/L001/L002} ${a/L001/L003} ${a/L001/L004} > ../merged_fastq/${BASENAME}
	done

for a in $(ls individual_fastq/*L001_R2_001.fastq.gz)
	do
		BASENAME=$(basename ${a})
		cat ${a} ${a/L001/L002} ${a/L001/L003} ${a/L001/L004} > ../merged_fastq/${BASENAME}
	done
	

## untrimmed fastqc -- typically skip this
#for a in $(find merged_fastq -type f \( -name "*.fastq.gz" \))
#	do
# 			echo Running fastqc on ${a}
#			fastqc -t 2 --extract -o fastqc_untrimmed_out ${a} 
#			echo ***Completed***
#	done

# trim and generate qc reports. RUN SCRIPT IN TRIMMED_FASTQ output directory
for a in $(ls merged_fastq/*R1_001.fastq.gz) 
	do
			echo ${a}
			echo ${a/R1/R2}

			/Users/nataliethornton/biotools/TrimGalore-0.6.5/trim_galore \
					--path_to_cutadapt /Users/nataliethornton/anaconda3/envs/wgs/bin/cutadapt \
					--quality 30 \
					--paired \
					--fastqc \
					"${a}" "${a/R1/R2}"
	done

#					--output_dir ~/seds/wgs/wgs_data/temp/temp_trimmed trimmed_fastq \
#					--fastqc_args "-o ~/seds/wgs/wgs_data/temp/temp_fastqc_trimmed/" \


