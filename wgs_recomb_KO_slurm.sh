#!/bin/bash -l
#SBATCH --partition=panda   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=CT_pckA_KO
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
IN_DIR=/athena/schnappingerlab/scratch/kar4019/wgs/23447Ehr/correct_23447Ehr/23447Ehr_N23242/caro/trimmed_fastq/
REF_FA=/athena/schnappingerlab/scratch/kar4019/wgs/23447Ehr/correct_23447Ehr/23447Ehr_N23242/caro/pckA_KO/pckA_KO_erdman_ref.fa
REF_GB=/athena/schnappingerlab/scratch/kar4019/Reference/Erdman/Erdman.gbk

Ref_strain=Erdman
Lab_strain=Erdman
PDIM_strain=Erdman

SAMPLE_LIST=/athena/schnappingerlab/scratch/kar4019/wgs/23447Ehr/correct_23447Ehr/23447Ehr_N23242/caro/sample_name_list_erdman
R1_FILE_ENDING=R1_001_val_1.fq.gz
R2_FILE_ENDING=R2_001_val_2.fq.gz

#spack load /3jymtx6 #samtools
#spack load bwa@0.7.17%gcc@8.2.0 arch=linux-centos7-sandybridge
#
### run alignment: 
#sh ~/wgs_scripts_dev/bash_scripts/align_and_clean_specific.sh \
#		-R "$REF_FA" -I "$IN_DIR" -N "$SLURM_CPUS_PER_TASK" -F "$SAMPLE_LIST" \
#		-1 "$R1_FILE_ENDING" -2 "$R2_FILE_ENDING"
#
## run SNP caller
sh ~/wgs_scripts_dev/bash_scripts/call_snv.sh -R "$REF_FA" 

spack unload /3jymtx6
spack unload bwa@0.7.17%gcc@8.2.0 arch=linux-centos7-sandybridge

# set everything up for running snp_csv_annotation.py script:
spack load blast-plus@2.9.0%gcc@8.2.0 arch=linux-centos7-sandybridge
export PATH="${python_path}:${PATH}"

echo "%%% PRINTING PATH"
echo $PATH

## run vcf parsing
sh ~/wgs_scripts_dev/bash_scripts/parse_vcf.sh -R "$REF_FA"

# DIRECTORY SETUP
mkdir snp_xlsx
mkdir blast

python ~/wgs_scripts_dev/python_analysis_scripts/snp_csv_annotation.py \
	   	-ref_strain recombinant \
		--H37Rv \
		-blast blast \
		-lab_strain "$Lab_strain" \
		-vcf vcf_parse \
		-out snp_xlsx

sh ~/wgs_scripts_dev/bash_scripts/pdim_check.sh \
		-S "$PDIM_strain" \
		-R $REF_GB \
		-C snp_xlsx

#make bedgraph files:
mkdir bedgraph

for f in bam/*dedup.bam
do
		echo "Writing .bedgraph for ${f}"
		BASE=$(basename ${f})

		bedtools genomecov -ibam "$f" -bg > bedgraph/"${BASE/bam/bedgraph}"
done

exit
