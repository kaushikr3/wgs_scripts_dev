#!/bin/bash

ref_genome='../../../Reference/H37RvCO/H37RvCO.gbk'

# mkdir cov
# for a in $(ls ./bam/*dedup.bam)
	# do
			# BASENAME=$(basename $a)
			# echo Running bedtools on $BASENAME
			# bedtools genomecov -ibam $a -d > cov/${BASENAME/_dedup.bam/.cov}
	# done

echo Running pdim_check
python ../../../Scripts/pdim_check.py -strain H37RvCO -stringency stringent -genome $ref_genome