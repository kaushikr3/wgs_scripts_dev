#!/bin/bash -l
#SBATCH --partition=panda   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=wt
#SBATCH --time=24:00:00   # HH/MM/SS
#SBATCH --mem=32G   # memory requested, units available: K,M,G,T
 
source ~/.bashrc
conda activate wgs 

python_path=$(dirname $(which python))

echo "This is job #:" $SLURM_JOB_ID >> slurm_output.log
echo "Running on node:" `hostname` >> slurm_output.log
echo "CPUS per task assigned:" "$SLURM_CPUS_PER_TASK" >> slurm_output.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> slurm_output.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> slurm_output.log


## VARIABLES TO SET RUN-WISE::
IN_DIR=/athena/schnappingerlab/scratch/kar4019/trimmed_fastq/
REF_FA=/athena/schnappingerlab/scratch/kar4019/wgs/recom_ref_DL_OFF
REF_GB=/athena/schnappingerlab/scratch/kar4019/Reference/H37RvCO/H37RvCO.gbk

Ref_strain=H37Rv

SAMPLE_LIST=/athena/schnappingerlab/scratch/kar4019/sample_name_list
R1_FILE_ENDING=R1_001_val_1.fq.gz
R2_FILE_ENDING=R2_001_val_2.fq.gz

#spack load -r python@3.7.0^gcc@6.3.0
spack load /3jymtx6 #samtools
spack load bwa@0.7.17%gcc@8.2.0 arch=linux-centos7-sandybridge
spack load bedtools2@2.28.0%gcc@6.3.0

## run alignment: 
sh ~/wgs_scripts_dev/bash_scripts/align_and_clean_specific.sh \
		-R "$REF_FA" -I "$IN_DIR" -N "$SLURM_CPUS_PER_TASK" -F "$SAMPLE_LIST" \
		-1 "$R1_FILE_ENDING" -2 "$R2_FILE_ENDING

##make bedgraph files:
mkdir bedgraph


for f in bam/*dedup.bam
do
		echo "Writing .bedgraph for ${f}"
		BASE=$(basename ${f})

		bedtools genomecov -ibam "$f" -bg > bedgraph/"${BASE/bam/bedgraph}"
done

exit
