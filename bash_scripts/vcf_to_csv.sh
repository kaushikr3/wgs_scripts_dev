#!/bin/bash

mkdir csv
mkdir blast

ref_genome='../../../Reference/H37RvCO/H37RvCO.fasta'
anno_genome='../../../Reference/AL123456/AL123456.gbk'
lab_ref_csv_dir='../../Lab_references/H37Rv_Jamie_North/csv/'

for a in $(ls ./vcf/*haploid.stringent.vcf)
	do
 			echo Running vcf_to_csv on $a
			python ../../../Scripts/vcf_to_csv.py -vcf $a -ref_genome $ref_genome -anno_genome $anno_genome -lab_ref_csv_dir $lab_ref_csv_dir -csv_main_dir_name csv -blast_main_dir_name blast
			echo ***Completed***
	done

	
for a in $(ls ./vcf/*diploid.lenient.vcf)
	do
 			echo Running vcf_to_csv on $a
			python ../../../Scripts/vcf_to_csv.py -vcf $a -ref_genome $ref_genome -anno_genome $anno_genome -lab_ref_csv_dir $lab_ref_csv_dir -csv_main_dir_name csv -blast_main_dir_name blast
			echo ***Completed***
	done
	
