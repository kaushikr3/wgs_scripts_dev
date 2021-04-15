#!/bin/bash -l
#SBATCH --partition=panda   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=snippy
#SBATCH --time=12:00:00   # HH/MM/SS
#SBATCH --mem=32G   # memory requested, units available: K,M,G,T
 
source ~/.bashrc
 
echo "This is job #:" $SLURM_JOB_ID >> slurm_output.log
echo "Running on node:" `hostname` >> slurm_output.log
echo "CPUS per task assigned:" "$SLURM_CPUS_PER_TASK" >> slurm_output.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> slurm_output.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> slurm_output.log
 

# DEPENDENCIES:
spack load miniconda3@4.6.14
conda activate snippy

# DIRECTORY SETUP
mkdir snippy

echo " ----  Running SNIPPY with reference: ${R} ----"

# SNIPPY ANALYSIS
for f in fastq/*R1_001_val_1.fq.gz
	do
		   	BASE=$(basename ${f}) 
			echo "$R"
			
			snippy --cpus "$SLURM_CPUS_PER_TASK" --cleanup --force \
					--outdir snippy -prefix "${BASE/_R1_001_val_1.fq.gz/}" \
					--reference "$R" \
					--R1 "$f" --R2 "${f/R1_001_val_1/R2_001_val_2}"

	done


echo "Left aligning and parsing snippy VCFs"
for f in snippy/*.filt.vcf
do
		BASE=$(basename ${f})
		vcfleftalign -r "$R" "$f" > "${f/.vcf/.leftalign.vcf}"
		
		~/biotools/gatk-4.2.0.0/gatk VariantsToTable -V "${f/.vcf/.leftalign.vcf}" \
				-F CHROM -F POS -F REF -F ALT -F QUAL -F AC -F AF -F QD \
				-GF AD -GF DP -GF GQ -GF GT -GF PL \
				-O vcf_parse/"${BASE/.vcf/_snippy.tsv}"

done


conda deactivate
spack unload miniconda3@4.6.14


