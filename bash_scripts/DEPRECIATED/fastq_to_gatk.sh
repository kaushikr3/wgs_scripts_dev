#!/bin/bash

#Run this script with sudo -E env "PATH=$PATH" ./file_name.sh

reference_genome=~/seds/ref_genomes/h37RvCO/h37RvCO4.fa
echo $PATH
mkdir bam
mkdir fastqc_untrimmed_out
mkdir sam
mkdir vcf
mkdir merged_fastq

for a in $(ls individual_fastq/*L001_R1_001.fastq.gz)
	do
		BASENAME=$(basename ${a})
		cat ${a} ${a/L001/L002} ${a/L001/L003} ${a/L001/L004} > merged_fastq/${BASENAME}
	done

for a in $(ls individual_fastq/*L001_R2_001.fastq.gz)
	do
		BASENAME=$(basename ${a})
		cat ${a} ${a/L001/L002} ${a/L001/L003} ${a/L001/L004} > merged_fastq/${BASENAME}
	done
	
	
#Change java version to 11 for fastqc
#update-alternatives --set java /usr/lib/jvm/java-11-openjdk-amd64/bin/java 

#for a in $(find './merged_fastq/*.fastq.gz')
for a in $(find merged_fastq -type f \( -name "*.fastq.gz" \))
	do
 			echo Running fastqc on ${a}
			fastqc -t 2 --extract -o fastqc_untrimmed_out ${a} 
			echo ***Completed***
	done



for a in $(ls ./merged_fastq/*R1_001.fastq.gz) 
	do

			# Use MC2155 for WGS of Msmeg
			# Use H37RvCO4 for WGS of H37Rv Mtb to get closest reference genome
			BASENAME=$(basename ${a})
			echo Aligning ${BASENAME} and ${BASENAME/R1/R2} with BWA MEM
			
			bwa mem -M -t 8 ${reference_genome} ${a} ${a/R1/R2} > sam/${BASENAME/R1_001.fastq.gz/bwa.aligned.sam}

			echo ${BASENAME/R1_001.fastq.gz/bwa.aligned.sam} is sorting with samtools

			# u is uncompressed BAM output,I think | is for chaining, -M is for marking shorter split hits as secondary-this is for picard compatibility
			samtools view -u sam/${BASENAME/R1_001.fastq.gz/bwa.aligned.sam} | samtools sort -o bam/${BASENAME/R1_001.fastq.gz/sorted.bam} -O bam -T temp -@ 8 -
			echo Wrote ${BASENAME/R1_001.fastq.gz/sorted.bam}
			rm sam/${BASENAME/R1_001.fastq.gz/bwa.aligned.sam}
			echo Removed ${BASENAME/R1_001.fastq.gz/bwa.aligned.sam}
			echo  
	done
	
		


jenv shell openjdk64-1.8.0.265

for a in ./bam/*sorted.bam
	do
	
			echo ${a} is running picard
			java -jar ~/biotools/picard/build/libs/picard.jar AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT I=${a} O=${a/sorted/RG} RGID=${a/_sorted.bam/} RGLB=${a/_sorted.bam/} RGPL=ILLUMINA RGPU=${a/_sorted.bam/} RGSM=${a/_sorted.bam/}
			java -jar ~/biotools/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT INPUT=${a/sorted/RG} OUTPUT=${a/sorted/dedup} METRICS_FILE=${a/sorted.bam/dedup.metrics.txt} USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
	done    
	


for a in ./bam/*dedup.bam 
	do
			echo ${a} is running samtools	 
			samtools index ${a}	 
			echo $a >> samtools_stats.log; samtools flagstat $a >> samtools_stats.log; echo %%%%%%%%%%%%%%% >> samtools_stats.log; echo  >> samtools_stats.log
	done
	
	
	
#update-alternatives --set java /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java	
#jenv shell openjdk64-1.8.0.265

for a in ./bam/*dedup.bam
	do
	BASENAME=$(basename $a)
	echo Running $BASENAME through GATK
	
	java -jar ~/biotools/GenomeAnalysisTK.jar -T HaplotypeCaller -stand_call_conf 30 -stand_emit_conf 10 -R $reference_genome -I $a --sample_ploidy 1 -o vcf/${BASENAME/dedup.bam/haploid.stringent.vcf}
	java -jar ~/biotools/GenomeAnalysisTK.jar -T HaplotypeCaller -stand_call_conf 3 -stand_emit_conf 3 -R $reference_genome -I $a --sample_ploidy 2 -o vcf/${BASENAME/dedup.bam/diploid.lenient.vcf}
	#java -jar ~/biotools/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar HaplotypeCaller -stand-call-conf 30 --sample-ploidy 1 -R $reference_genome -I $a -O vcf/${BASENAME/dedup.bam/haploid.stringent.NEW.vcf}
	#java -jar ~/biotools/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar HaplotypeCaller -stand-call-conf 3 --sample-ploidy 2 -R $reference_genome -I $a -O vcf/${BASENAME/dedup.bam/diploid.lenient.NEW.vcf}

	done
	

	

