### WGS scripts dev'd 2020-2022

General workflow:
1. Trim and run any necessary QC on concerning samples (trimming code in trim_and_qc_slurm.sh)
    - Set filelist name to path to file holding filepaths for all R1 raw fastq files to be trimmed.
    - Alter out_dir name ('trimmed_fastq') if desired
3. Generate any recombinant reference genomes requried (currently peformed elsewhere, on my local setup)
4. Run code as shown in sample wgs_slurm.sh file. Key steps:
    - Alter SLURM resource calls as needed
    - Set reference genome path for 'REF_FA' and 'REF_GB'
    - Set file source directory (by default will be 'trimmed_fastq')
    - Set ref strain based on experimental/species strains. Set={'H37Rv', 'Erdman', 'HN878', 'BCG', 'Msmeg', 'AL123456', 'recombinant'}
    - Set lab strain based on lab strain used, if lab SNP reference exists. Else use None. Set={'North', 'UMass', 'Erdman', 'HN878', 'BCG', None}
    - Set ref strain for PDIM reference. Set={'H37RvCO', 'HN878', 'Erdman', 'BCG', 'AL123456', 'NC_000962'}
    - Generate a file holding a list of R1 fastq file names and assign its name to 'FILELIST' variable. This list will just contain filenames, the directory path is passed in as 'IN_DIR'
    - Script will run through:
        - align_and_clean.sh
            a. Runs BWA alignment
            b. Runs Picard tools CleanSam
            c. Runs Picard tools AddOrReplaceReadGroups
            d. sorts .bam files
            e. Runs Picard tools MarkDuplicates
            f. indexes .bam files
            g. generates flagstat report file on .bam files
        - call_snv.sh
            a. Runs GATK HaplotypeCaller (quasi-Stringent; done by running on --sample-ploidy 1)
            b. Runs GATK HaplotypeCaller (quasi-Lenient; done by running on --sample-ploidy 2; which relaxes SNP thresholding)
        - parse_vcf.sh
            Parses GATK file to .tsv
        - snp_csv_annotation.sh
            Generates an .xlsx file for each sample consisting of an annotated stringent and lenient SNP sheet
        - pdim_check.sh
            a. Generates a coverage file from each sorted, deduplicates .bam file
            b. Generates an .xlsx file holding PDIM status for each sample based on SNPs and coveragec
        - call_structural_variants.sh
            Runs pilon to get structural analysis of sample
5. After running, you can generate .bedgraph files from sorted, deduplicated .bam files for visualization of read coverage as needed using 'gen_bedgraph_slurm.sh'


NOTES: 
Do not run PDIM check using recombinant reference alignment files -- this script is reference coordinate dependant

          
