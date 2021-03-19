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
        echo '   -ref   path to reference_genome.fa'
        echo '   -in   path to directory holding UNMERGED fastq'
        exit 1
}

if [[ ${#} -eq 0 ]]; then
   usage
fi

# Define list of arguments expected in the input
optstring="ref:in:"

while getopts ${optstring} arg; do
  case "${arg}" in
    ref) REF=${OPTARG} ;;
    in) UNMERGED_DIR=${OPTARG} ;;

    ?)
      echo "Invalid option: -${OPTARG}."
      echo
      usage
      ;;
  esac
done

echo " ----  Processing fastq.gz files from ${UNMERGED_DIR}  ----" >> reports/"${UNMERGED_DIR}.log"
echo " ----  Aligning files against ref: ${REF} ----" >> reports/"${UNMERGED_DIR}.log"


# DIRECTORY SETUP
mkdir fastq
mkdir fastq/merged_untrimmed

mkdir reports
mkdir reports/fastqc_untrimmed_out
mkdir reports/fastqc_trimmed_out

mkdir bam
mkdir gatk
mkdir freebayes


# MERGE LANE FASTQS
for f in ${UNMERGED_DIR}/*L001_R1_001.fastq.gz
do
  [[ -e "$f" ]] || break  # handle the case of no *L001_R1_001.fastq.gz files

  BASE="${f/L001_R1_001.fastq.gz/}"
  FILES_R1=(${BASE}*R1_001.fastq.gz)
  cat "${FILES_R1[@]}" > fastq/merged_untrimmed/"${BASE}R1_001.fastq.gz"
  echo "Merged "${FILES_R1}" into ${BASE}_R1_001.fastq.gz" >> reports/"${UNMERGED_DIR}.log"

  FILES_R2=(${BASE}*R2_001.fastq.gz)
  cat "${FILES_R2[@]}" > fastq/merged_untrimmed/"${BASE}R2_001.fastq.gz"
  echo "Merged "${FILES_R2}" into ${BASE}_R2_001.fastq.gz" >> reports/"${UNMERGED_DIR}.log"

done


# UNTRIMMED FASTQC
for f in fastq/merged_untrimmed/*_001.fastq.gz

	do
 			echo Running fastqc on "$f" >> reports/"${UNMERGED_DIR}.log"
			fastqc -t 2 --extract -o reports/fastqc_untrimmed_out "$f"
	done


# TRIM AND GENERATE TRIMMED FASTQC
for f in fastq/merged_untrimmed/*R1_001.fastq.gz
	do
			echo "Running TrimGalore on:  ${a}  ${a/R1/R2}" >> reports/"${UNMERGED_DIR}.log"
			/Users/nataliethornton/biotools/TrimGalore-0.6.5/trim_galore \
					--path_to_cutadapt /Users/nataliethornton/anaconda3/envs/wgs/bin/cutadapt \
					--quality 30 \
					--paired \
					--fastqc \
					--output_dir fastq \
					--fastqc_args "-t 2 --extract -o reports/fastqc_trimmed_out" \
					"${f}" "${f/R1/R2}"
	done


# ALIGNMENT:
for f in fastq/*R1_001_val_1.fq.gz
	do
			# Use MC2155 for WGS of Msmeg
			# Use H37RvCO4 for WGS of H37Rv Mtb to get closest reference genome
			echo "Aligning ${f} and ${f/R1_001_val_1/R2_001_val_2} with BWA MEM" >> reports/"${UNMERGED_DIR}.log"

			bwa mem -M -t 8 "$REF" "$f" "${f/R1_001_val_1/R2_001_val_2}" | \
			    samtools view -Shb - > bam/"${f/R1_001_val_1.fq.gz/unsorted.bam}"

	done


# PICARD TOOLS
for f in bam/*unsorted.bam
	do

			echo "Cleaning, Running PICARD Tools, and sorting ${f} " >> reports/"${UNMERGED_DIR}.log"

			# clean unsorted bamfile:
			java -jar $picard_path CleanSam \
			    INPUT=bam/"$f" \
			    OUTPUT=bam/"${f/unsorted.bam/unsorted.cleaned.bam}"

      # sort bam
      samtools sort -o bam/"${f/unsorted/sorted}" -O bam -T temp -@ 8 bam/"${f/unsorted.bam/unsorted.cleaned.bam}"

      # removing the add or replace read groups steps here
#			java -jar ~/biotools/picard/build/libs/picard.jar AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT \
#					I="$f" O="${f/sorted/RG}" RGID="${f/_sorted.bam/}" RGLB="${a/_sorted.bam/}" RGPL=ILLUMINA \
#					RGPU="${f/_sorted.bam/}" RGSM="${f/_sorted.bam/}" USE_JDK_DEFLATER=true USE_JDK_INFLATER=true

			java -jar $picard_path MarkDuplicates \
          VALIDATION_STRINGENCY=STRICT \
					INPUT=bam/"${f/unsorted/sorted}" \
					OUTPUT=bam/"${f/unsorted/strict.dedup}" \
					METRICS_FILE=reports/"${f/unsorted.bam/strict.dedup.metrics.txt}" \
			    ASSUME_SORTED=true

			java -jar ~/biotools/picard/build/libs/picard.jar MarkDuplicates \
					VALIDATION_STRINGENCY=LENIENT \
					INPUT=bam/"${f/unsorted/sorted}" \
					OUTPUT="${f/unsorted/lenient.dedup}" \
					METRICS_FILE=reports/"${f/sorted.bam/lenient.dedup.metrics.txt}" \
          ASSUME_SORTED=true
#					REMOVE_DUPLICATES=TRUE \
#					USE_JDK_DEFLATER=true \
#					USE_JDK_INFLATER=true

	done


# INDEX BAMS
for f in bam/*dedup.bam
	do
			echo "Indexing deduplicated bam: ${f}" >> reports/"${UNMERGED_DIR}.log"
			samtools index "$f"
			echo "$f" >> reports/samtools_stats.log
			samtools flagstat "$f" >> reports/samtools_stats.log
			echo %%%%%%%%%%%%%%% >> reports/samtools_stats.log
			echo  >> reports/samtools_stats.log
	done


# GATK ANALYSIS
for f in bam/*dedup.bam
	do
	echo "Running GATK on ${f}" >> reports/"${UNMERGED_DIR}.log"

	java -jar ~/biotools/gatk-3.5/GenomeAnalysisTK.jar -T HaplotypeCaller \
	  -R "$REF" -I "$f" \
	  --sample_ploidy 1 \
	  -O gatk/"${f/dedup.bam/gatk.haploid.vcf}"

  java -jar $gatk_path HaplotypeCaller \
    -R "$REF" -I "$f" \
    -O gatk/"${f/dedup.bam/gatk.diploid.vcf}"

	done


# FREEBAYES
for f in bam/*.dedup.bam
do
  echo "Running freebayes on ${f}" >> reports/"${UNMERGED_DIR}.log"
  freebayes -f "$REF" --ploidy 1 bam/"$f" > freebayes/"${f/dedup.bam/freebayes.vcf}"



## SNIPPY
#  $snippy_path --cpus $num_procs --outdir $out_dir/$root-snippy --prefix $root --cleanup --ref $ref --R1 $fq_1 --R2 $fq_2
#  print SH "mv $out_dir/$root-snippy/$root.filt.vcf $out_dir/$root.$aligner.$caller.vcf
#  rm -r $out_dir/$root-snippy