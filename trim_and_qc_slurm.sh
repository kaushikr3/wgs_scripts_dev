#!/bin/bash -l
#SBATCH --partition=panda   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=yuan_trim
#SBATCH --time=12:00:00   # HH/MM/SS
#SBATCH --mem=32G   # memory requested, units available: K,M,G,T
 
source ~/.bashrc
conda activate trim 


echo "This is job #:" $SLURM_JOB_ID >> slurm_output.log
echo "Running on node:" `hostname` >> slurm_output.log
echo "CPUS per task assigned:" "$SLURM_CPUS_PER_TASK" >> slurm_output.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> slurm_output.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> slurm_output.log

spack load /v46thbu #fastqc
spack load /4mpgaad #py-cutadapt2.10 

filelist=raw_fastq_list
out_dir=trimmed_fastq

mkdir "$out_dir"

while read f 
	do
	/home/kar4019/biotools/TrimGalore-0.6.10/trim_galore --paired --output_dir "$out_dir" "${f}" "${f/_R1_001.fastq.gz/_R2_001.fastq.gz}"

done< "$filelist"


mkdir qc_reports

for f in trimmed_fastq/*val_1.fq.gz 
do
		fastqc -o qc_reports/ "${f}" "${f/1_val_1/2_val_2}"

done

exit
