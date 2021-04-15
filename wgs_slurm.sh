#!/bin/bash -l
#SBATCH --partition=panda   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=wgs_align_and_call_snv
#SBATCH --time=12:00:00   # HH/MM/SS
#SBATCH --mem=32G   # memory requested, units available: K,M,G,T
 
source ~/.bashrc
 
echo "This is job #:" $SLURM_JOB_ID >> slurm_output.log
echo "Running on node:" `hostname` >> slurm_output.log
echo "CPUS per task assigned:" "$SLURM_CPUS_PER_TASK" >> slurm_output.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> slurm_output.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> slurm_output.log
 

source ~/wgs/wgs_scripts_dev/run_full_pipeline.sh \
		-R ~/wgs/Reference/{PATH/TO/REF.fasta} \
		-G ~/wgs/Reference/{PATH/TO/REF.gbk} \
		-I /athena/schnappingerlab/scratch/nat4004/{PATH/TO/UNMERGED/FASTQ}/ \
		-N "$SLURM_CPUS_PER_TASK"

exit

