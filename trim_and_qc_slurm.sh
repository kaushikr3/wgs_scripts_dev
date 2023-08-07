#!/bin/bash -l
#SBATCH --partition=panda   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=trim3
#SBATCH --time=12:00:00   # HH/MM/SS
#SBATCH --mem=32G   # memory requested, units available: K,M,G,T
 
source ~/.bashrc
 
echo "This is job #:" $SLURM_JOB_ID >> slurm_output.log
echo "Running on node:" `hostname` >> slurm_output.log
echo "CPUS per task assigned:" "$SLURM_CPUS_PER_TASK" >> slurm_output.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> slurm_output.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> slurm_output.log

spack load fastqc
spack load -r py-cutadapt@1.13
spack load trimgalore

filelist=raw_fastq_list
out_dir=trimmed_fastq

mkdir "$out_dir"

while read f 
	do
			trim_galore --paired --output_dir "$out_dir" \
					raw_fastq/"${f}" raw_fastq/"${f/_1\.fastq/_2\.fastq}"
	done< "$filelist"


mkdir qc_reports

for f in trimmed_fastq/*val_1.fq.gz 
do
		fastqc -o qc_reports/ "${f}" "${f/1_val_1/2_val_2}"

done

exit
