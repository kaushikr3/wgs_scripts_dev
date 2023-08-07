#!/bin/bash

# DEPENDENCIES:
# fastqc
# trimgalore
# cutadapt
# bwa
# samtools
# picard
# gatk
# freebayes


# READ IN ARGUMENTS
# usage:
function usage {
        echo "Usage: $(basename $0) [-ref] [-in]" 2>&1
		        echo '   -R   path to reference_genome.fa'
				        echo '   -I   path to directory holding UNMERGED fastq'
						        exit 1
						}

						if [[ ${#} -eq 0 ]]; then
								   usage
						   fi

						   # Define list of arguments expected in the input
						   optstring="R:I:"

						   while getopts ${optstring} arg; do
								     case "${arg}" in
											     R) REF="${OPTARG}" ;;
												     I) UNMERGED_DIR="${OPTARG}" ;;

													     ?)
																       echo "Invalid option: -${OPTARG}."
																	         echo
																			       usage
																				         ;;
																						   esac
																				   done

																				   # DIRECTORY SETUP
																				   mkdir fastq
																				   mkdir fastq/merged_untrimmed

																				   mkdir reports
																				   mkdir reports/fastqc_untrimmed_out
																				   mkdir reports/fastqc_trimmed_out

																				   mkdir bam
																				   mkdir gatk
																				   mkdir freebayes


																				   echo " ----  Processing fastq.gz files from ${UNMERGED_DIR}  ----"
																				   echo " ----  Aligning files against ref: ${REF} ----"


																				   # MERGE LANE FASTQS
																				   for f in ${UNMERGED_DIR}/*L001_R1_001.fastq.gz
																				   do 
																						   		#[[ -e "$f" ]] || break  # handle the case of no *L001_R1_001.fastq.gz files
																									
																										echo "$f"
																										  		BASE="${f/L001_R1_001.fastq.gz/}"
																														OUTNAME=$(basename ${BASE})

																														  		FILES_R1=(${BASE}*R1_001.fastq.gz)
																																  		FILES_R2=(${BASE}*R2_001.fastq.gz)
																																				
																																		  		cat "${FILES_R1[@]}" > fastq/merged_untrimmed/"${OUTNAME}R1_001.fastq.gz"
																																				  		echo "Merged "${FILES_R1}" into ${OUTNAME}_R1_001.fastq.gz"

																																						  		cat "${FILES_R2[@]}" > fastq/merged_untrimmed/"${OUTNAME}R2_001.fastq.gz"
																																								  		echo "Merged "${FILES_R2}" into ${OUTNAME}_R2_001.fastq.gz"

																																								done


																																								# UNTRIMMED FASTQC
																																								for f in fastq/merged_untrimmed/*_001.fastq.gz

																																											do
																																													 			echo Running fastqc on "$f"
																																																			fastqc -t 2 --extract -o reports/fastqc_untrimmed_out "$f"
																																																				done


																																																				# TRIM AND GENERATE TRIMMED FASTQC
																																																				for f in fastq/merged_untrimmed/*R1_001.fastq.gz
																																																							do
																																																												echo "Running TrimGalore on:  ${f}  ${f/R1/R2}"
																																																															trim_galore --quality 30 --paired --output_dir fastq \
																																																																						--fastqc --fastqc_args "-t 2 --extract -o reports/fastqc_trimmed_out" \
																																																																											"${f}" "${f/R1/R2}"
																																																																			#	--path_to_cutadapt /Users/nataliethornton/anaconda3/envs/wgs/bin/cutadapt \
																																																																						done


																																																																			# ALIGNMENT:
																																																																			for f in fastq/*R1_001_val_1.fq.gz
																																																																						do
																																																																											# Use MC2155 for WGS of Msmeg
																																																																														# Use H37RvCO4 for WGS of H37Rv Mtb to get closest reference genome
																																																																																	echo "Aligning ${f} and ${f/R1_001_val_1/R2_001_val_2} with BWA MEM" 

																																																																																				OUTNAME=$(basename ${f})

																																																																																							bwa mem -M -t 8 "$REF" "$f" "${f/R1_001_val_1/R2_001_val_2}" | \
																																																																																												    samtools view -Shb - > bam/"${OUTNAME/R1_001_val_1.fq.gz/unsorted.bam}"

																																																																																								done


																																																																																								# PICARD TOOLS
																																																																																								for f in bam/*unsorted.bam
																																																																																											do
																																																																																																echo "Cleaning and sorting ${f} " 
																																																																																																			BASE=$(basename ${f})
																																																																																																						
																																																																																																						# clean unsorted bamfile:
																																																																																																									picard CleanSam \
																																																																																																																INPUT="$f" \
																																																																																																																					OUTPUT="${f/unsorted.bam/unsorted.cleaned.bam}"
																																																																																																										
																																																																																																												# sort bam 
																																																																																																															samtools sort -o "${f/unsorted/sorted}" -O bam -T temp -@ 8 "${f/unsorted.bam/unsorted.cleaned.bam}"
																																																																																																																	   
																																																																																																																		echo "Running MarkDuplicates ${f}" 
																																																																																																																					picard MarkDuplicates \
																																																																																																																												VALIDATION_STRINGENCY=STRICT \
																																																																																																																																	INPUT="${f/unsorted/sorted}" \
																																																																																																																																						OUTPUT="${f/unsorted/strict.dedup}" \
																																																																																																																																											METRICS_FILE=reports/"${BASE/unsorted.bam/strict.dedup.metrics.txt}" \
																																																																																																																																																ASSUME_SORTED=true
																																																																																																																						
																																																																																																																						    	  # removing the add or replace read groups steps here
																																																																																																																								  					java -jar ~/biotools/picard/build/libs/picard.jar AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT \
																																																																																																																																						I="$f" O="${f/sorted/RG}" RGID="${f/_sorted.bam/}" RGLB="${a/_sorted.bam/}" RGPL=ILLUMINA \
																																																																																																																																													RGPU="${f/_sorted.bam/}" RGSM="${f/_sorted.bam/}" USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
																																																																																																																														
																																																																																																																																		#java -jar $picard_path MarkDuplicates \
																																																																																																																																									#java -jar ~/biotools/picard/build/libs/picard.jar MarkDuplicates \
																																																																																																																																										
																																																																																																																																												picard MarkDuplicates \ 
																																																																																																																																							VALIDATION_STRINGENCY=LENIENT \
																																																																																																																																														INPUT="${f/unsorted/sorted}" \
																																																																																																																																																			OUTPUT="${f/unsorted/lenient.dedup}" \
																																																																																																																																																								METRICS_FILE=reports/"${BASE/sorted.bam/lenient.dedup.metrics.txt}" \
																																																																																																																																																													ASSUME_SORTED=true
																																																																																																																																								#				REMOVE_DUPLICATES=TRUE \
																																																																																																																																											#				USE_JDK_DEFLATER=true \
																																																																																																																																												#				USE_JDK_INFLATER=true

																																																																																																																																											echo "Indexing ${f}" 
																																																																																																																																														samtools index "$f"
																																																																																																																																															
																																																																																																																																															done


																																																																																																																																															# EVAL BAMS
																																																																																																																																															for f in bam/*dedup.bam
																																																																																																																																																		do
																																																																																																																																																							echo "BAM  ${f}" >> reports/samtools_stats.log
																																																																																																																																																										samtools flagstat $a >> reports/samtools_stats.log
																																																																																																																																																													echo %%%%%%%%%%%%%%% >> reports/samtools_stats.log
																																																																																																																																																																echo  >> reports/samtools_stats.log
																																																																																																																																																																	done


																																																																																																																																																																	# GATK ANALYSIS
																																																																																																																																																																	for f in bam/*dedup.bam
																																																																																																																																																																				do
																																																																																																																																																																									echo "Running GATK on ${f}" 
																																																																																																																																																																												BASE=$(basename ${f})
																																																																																																																																																																															
																																																																																																																																																																															~/biotools/gatk-4.2.0.0/gatk HaplotypeCaller \
																																																																																																																																																																																						-R "$REF" -I "$f" \
																																																																																																																																																																																										   	--sample_ploidy 1 \
																																																																																																																																																																																															   	-O gatk/"${BASE/dedup.bam/gatk.haploid.vcf}"
																																																																																																																																																																																	   
																																																																																																																																																																																		~/biotools/gatk-4.2.0.0/gatk HaplotypeCaller \
																																																																																																																																																																																								   	-R "$REF" -I "$f" \
																																																																																																																																																																																													   	-O gatk/"${BASE/dedup.bam/gatk.diploid.vcf}"

																																																																																																																																																																																			done


																																																																																																																																																																																			# FREEBAYES

																																																																																																																																																																																			for f in bam/*.dedup.bam
																																																																																																																																																																																			do
																																																																																																																																																																																						   	echo "Running freebayes on ${f}" 
																																																																																																																																																																																								   	BASE=$(basename ${f})

																																																																																																																																																																																										   	~/biotools/freebayes-1.3.4-linux-static-AMD64 -f "$REF" --ploidy 1 "$f" > freebayes/"${BASE/dedup.bam/freebayes.vcf}"
																																																																																																																																																																																									done


																																																																																																																																																																																									## SNIPPY
																																																																																																																																																																																									#  $snippy_path --cpus $num_procs --outdir $out_dir/$root-snippy --prefix $root --cleanup --ref $ref --R1 $fq_1 --R2 $fq_2
																																																																																																																																																																																									#  print SH "mv $out_dir/$root-snippy/$root.filt.vcf $out_dir/$root.$aligner.$caller.vcf
																																																																																																																																																																																									#  rm -r $out_dir/$root-snippy


