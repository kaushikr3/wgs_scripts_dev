#!/bin/bash -l
#SBATCH --partition=panda   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=wt
#SBATCH --time=24:00:00   # HH/MM/SS
#SBATCH --mem=32G   # memory requested, units available: K,M,G,T
 
source ~/.bashrc
conda activate snippy

python_path=$(dirname $(which python))

echo "This is job #:" $SLURM_JOB_ID >> slurm_output.log
echo "Running on node:" `hostname` >> slurm_output.log
echo "CPUS per task assigned:" "$SLURM_CPUS_PER_TASK" >> slurm_output.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> slurm_output.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> slurm_output.log


## VARIABLES TO SET RUN-WISE::
IN_DIR=/athena/schnappingerlab/scratch/nat4004/AN00012283/trimmed_fastq
REF_FA=~/wgs/Reference/H37RvCO/H37RvCO.fasta
REF_GB=~/wgs/Reference/H37RvCO/H37RvCO.gbk

Ref_strain=H37Rv
Lab_strain=UMass
PDIM_strain=H37RvCO

FILELIST=filelist
BAMLIST=bamlist

spack load -r python@3.7.0^gcc@6.3.0
spack load samtools@1.9%gcc@6.3.0
spack load bwa@0.7.15%gcc@6.3.0

## run alignment: 
#sh ~/wgs/wgs_scripts_dev/bash_scripts/align_and_clean.sh -R "$REF_FA" -I "$IN_DIR" -N "$SLURM_CPUS_PER_TASK"
sh ~/wgs/wgs_scripts_dev/bash_scripts/align_and_clean_specific.sh \
		-R "$REF_FA" -I "$IN_DIR" -N "$SLURM_CPUS_PER_TASK" -F "$FILELIST"

## run SNP caller
#sh ~/wgs/wgs_scripts_dev/bash_scripts/call_snv.sh -R "$REF_FA" 
sh ~/wgs/wgs_scripts_dev/bash_scripts/call_snv_specific.sh -R "$REF_FA" -F "$BAMLIST"

spack unload samtools@1.9%gcc@6.3.0
spack unload bwa@0.7.15%gcc@6.3.0

## set everything up for running snp_csv_annotation.py script:
spack load blast-plus@2.10.0
export PATH="${python_path}:${PATH}"

echo "%%% PRINTING PATH"
echo $PATH

## run vcf parsing
sh ~/wgs/wgs_scripts_dev/bash_scripts/parse_vcf.sh -R "$REF_FA"


## DIRECTORY SETUP
mkdir snp_xlsx
mkdir blast

python ~/wgs/wgs_scripts_dev/python_analysis_scripts/snp_csv_annotation.py \
	   	-ref_strain "$Ref_strain" \
		--H37Rv \
		-blast blast \
		-lab_strain "$Lab_strain" \
		-vcf vcf_parse \
		-out snp_xlsx

sh ~/wgs/wgs_scripts_dev/bash_scripts/pdim_check.sh \
		-S "$PDIM_strain" \
		-R $REF_GB \
		-C snp_xlsx

## run structural variant caller
sh ~/wgs/wgs_scripts_dev/bash_scripts/call_structural_variants.sh -R "$REF_FA" -N "$SLURM_CPUS_PER_TASK"

exit

