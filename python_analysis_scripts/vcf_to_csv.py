#!/usr/bin/env python
# coding: utf-8
from __future__ import annotations

import argparse
import os
import sys
import glob
from collections import defaultdict

import numpy as np
import pandas as pd
import pyranges as pr
from Bio import SearchIO, SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main():
    my_parser = argparse.ArgumentParser(prog='vcf_to_csv')
    my_parser.add_argument('-vcf', required=True, dest='vcf_path', type=str)
    my_parser.add_argument('-ref_genome', required=True,
                           dest='ref_genome_path', type=str)
    my_parser.add_argument('-anno_genome', required=True,
                           dest='anno_genome_path', type=str)
    my_parser.add_argument('-lab_ref_csv_dir', required=False,
                           dest='lab_ref_csv_dir', type=str, default=None)
    my_parser.add_argument('-csv_main_dir_name', type=str, default='csv')
    my_parser.add_argument('-blast_main_dir_name', type=str, default='blast')
    args = my_parser.parse_args()
    
    cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "sample"]
    vcf = pd.read_csv(args.vcf_path, sep='\t', comment='#', names=cols)
    #vcf = pd.read_csv(args.vcf_path, sep='\t', header=28)
    # print(os.getcwd())
    anno_genome = SeqIO.read(args.anno_genome_path, format="genbank")
    print(anno_genome.name)

    ref_genome = SeqIO.read(args.ref_genome_path, format="fasta")
    reference_path = os.path.split(os.path.split(args.ref_genome_path)[0])[0]

    # reference_path = '../../../../../Reference'
    #print(" !!! TEST: !!!!")
    print(reference_path)

    selected_features = pd.read_csv(f"~/wgs/metadata_wgs/{anno_genome.name}_selected_features.csv")
    selected_features_pr = pr.PyRanges(selected_features)
    if args.lab_ref_csv_dir:
        print("%%%%% Lab strains pulled in %%%%%")
        # print(args.lab_ref_csv_dir+'*.haploid.stringent.annotated.csv')
        stringent_lab_ref = pd.read_csv(
            glob.glob(os.path.join(args.lab_ref_csv_dir,'*haploid.stringent.annotated.csv'))[0])
        lenient_lab_ref = pd.read_csv(
            glob.glob(os.path.join(args.lab_ref_csv_dir,'*diploid.lenient.annotated.csv'))[0])
    else:
        stringent_lab_ref = None
        lenient_lab_ref = None
    # print(selected_features['Feature_type'].unique())
    basename = os.path.basename(args.vcf_path)
    csv_name = basename.replace('.vcf', '.annotated.csv')

    # MADE SOME CHANGES HERE; SWITCHED EVERYTHING TO OCCUR IN TEMP BLAST FILES;
    blast_dir_name = args.blast_main_dir_name
    #blast_dir_name = f"{args.blast_main_dir_name}/{''.join(basename.split('.')[:-1])}/"
    #os.system(f"mkdir {blast_dir_name}")

    assert sys.version_info >= (3, 6)	
    output_df = vcf_to_csv(vcf, ref_genome, anno_genome,
                           selected_features_pr, stringent_lab_ref, lenient_lab_ref, blast_dir_name,
                           os.path.split(args.anno_genome_path)[0],
                           reference_path=reference_path)
    output_df['Alt_freq'] = output_df['Alt_freq'].fillna(0)
    output_df['Alt_freq'] = output_df['Alt_freq'].replace('NA', 0)
    # print(output_df['Alt_freq'].unique())
    output_df = output_df.sort_values(by='Alt_freq', ascending=False)
    output_df.to_csv(f"{args.csv_main_dir_name}/{csv_name}", index=False)
    return None


def vcf_to_csv(vcf: pd.DataFrame, ref_genome: SeqRecord, anno_genome: SeqRecord, 
				selected_features_pr: pr.PyRanges, stringent_lab_ref: pd.DataFrame, lenient_lab_ref: pd.DataFrame, 
                blast_dir_name: str, annotation_path: str, reference_path="~/wgs/Reference") -> pd.DataFrame:
    """Takes in a vcf file from GATK, adds columns with annotations of the mutations and other details"""
    # Added check to make sure the genotype fields is correct since was having strange error where some rows would only be GT:GQ:PL.
    # Decided to delete these rows and print warning
    #print("!!!! TEST BELOW ")
    #print(vcf.columns)
    
    n_bad_rows = (vcf["FORMAT"]!="GT:AD:DP:GQ:PL").sum()
    if n_bad_rows > 0:
        print(f"!!!!Deleting {n_bad_rows} rows due to incomplete genotype information!!!!")
        vcf = vcf[vcf["FORMAT"]=="GT:AD:DP:GQ:PL"]
    output_df = pd.DataFrame(columns=['Strain', 'Ref_genome'])
    output_df[['Ref_pos', 'Ref', 'Alt']] = vcf[['POS', 'REF', 'ALT']]
    output_df['Strain'] = vcf.columns[9].split('/')[-1]
    output_df['Ref_genome'] = ref_genome.name
    output_df['#Ref_reads'] = vcf.iloc[:, 9].str.split(':|,', expand=True)[
        1].astype(int)
    output_df['#Alt_reads'] = vcf.iloc[:, 9].str.split(':|,', expand=True)[
        2].astype(int)
    output_df['Coverage'] = output_df['#Ref_reads'] + output_df['#Alt_reads']
    output_df['Alt_freq'] = output_df['#Alt_reads']/output_df['Coverage']
    output_df['GQ'] = vcf.iloc[:, 9].str.split(':|,', expand=True)[
        4].astype(int)
    output_df['Anno_genome'] = anno_genome.name

    #print(type(anno_genome))
    #print(type(ref_genome))
    #print(type(blast_dir_name))

	# get reference blast db directory:
    if os.path.isfile(os.path.join(annotation_path, anno_genome.name + ".fasta")):
        reference_location = os.path.join(annotation_path, anno_genome.name + ".fasta")
    
    elif os.path.isfile(os.path.join(annotation_path, anno_genome.name + ".fa")):
        reference_location = os.path.join(annotation_path, anno_genome.name + ".fa")
    
    elif os.path.isfile(os.path.join(annotation_path, os.path.split(annotation_path)[1] + ".fasta")):
        reference_location = os.path.join(annotation_path, os.path.split(annotation_path)[1] + ".fasta")
    
    elif os.path.isfile(os.path.join(annotation_path, os.path.split(annotation_path)[1] + ".fa")):
        reference_location = os.path.join(annotation_path, os.path.split(annotation_path)[1] + ".fa")
    
    else:
        reference_location = os.path.join(annotation_path, os.path.split(annotation_path)[1] + ".fasta")
        print(reference_location)
        assert os.path.isfile(reference_location) # blast reference db files not found

    output_df[['#HSPs', 'Percent_identity', 'Anno_position']] = vcf.apply(
        blast_ref, axis=1, args=(ref_genome, anno_genome, blast_dir_name, reference_location))
    output_df[['Location', 'Feature_type', 'Feature_name', 'Distance']] = output_df.apply(
        extract_features, axis=1, args=(selected_features_pr,))
    cols_to_explode = ['Feature_type', 'Feature_name', 'Distance']
    output_df = output_df.set_index([c for c in output_df.columns if c not in cols_to_explode]).apply(
        pd.Series.explode).reset_index()
    output_df[['AA_change', 'Mutation_type']
              ] = output_df.apply(annotate_aa, axis=1, args=(selected_features_pr, anno_genome))
    output_df['Distance'] = output_df['Distance'].astype(str)

    # set defaults for if metadata gene_detail sheets do not exist
    # gene_details = False
    # cols_to_group = ['Distance', 'Feature_type', 'Feature_name', 'AA_change', 'Mutation_type']

    # get metadata, if exists
    # if os.path.exists(os.path.join('~/seds/wgs/wgs_scripts/metadata/', anno_genome.name + '_gene_details.csv')):
    gene_details = pd.read_csv(
        f"~/wgs/metadata_wgs/{anno_genome.name}_gene_details.csv")
    output_df = pd.merge(output_df, gene_details, left_on=[
                         'Feature_name'], right_on=['ORF'], how='left').drop(columns='ORF')
    cols_to_group = ['Distance', 'Feature_type', 'Feature_name',
                     'Feature_gene_name', 'Feature_description', 'AA_change', 'Mutation_type']

    # else:
    #     print("WARNING: no gene details file found")

    output_df = output_df.fillna('NA')
    # group duplicates
    output_df = output_df.groupby([c for c in output_df.columns if c not in cols_to_group])[
        cols_to_group].agg(list).reset_index()
    output_df[cols_to_group] = output_df[cols_to_group].applymap(';'.join)
    # add lab_reference_info
    if (stringent_lab_ref is not None) and (lenient_lab_ref is not None):
        assert stringent_lab_ref['Strain'].values[0] == lenient_lab_ref['Strain'].values[0]
        # print(lab_ref)
        output_df['Lab_reference'] = stringent_lab_ref['Strain'].values[0]
        output_df['Present_in_lab_reference_stringent'] = output_df.apply(
            add_lab_ref, axis=1, args=(stringent_lab_ref,))
        output_df['Present_in_lab_reference_lenient'] = output_df.apply(
            add_lab_ref, axis=1, args=(lenient_lab_ref,))
    return output_df


def blast_ref(row: pd.Series, ref_genome: SeqRecord, anno_genome: SeqRecord, blast_dir_name: str, 
				reference_location: str, query_len=100) -> pd.Series:
    """For each row, uses the reference position to create a query sequence of length query length. This is then blast-ed against the annotated genome. The alignment is verified and details are returned as a Series"""
    # Note that reference genome is zero-indexed and ref position is 1-indexed
    ref_position = row['POS']
    ref = row['REF']
	
    #reference_location = os.path.join(reference_path, anno_genome.name, anno_genome.name + ".fasta")
    #if os.path.isfile(os.path.join(annotation_path, anno_genome.name + ".fasta")):
    #    reference_location = os.path.join(annotation_path, anno_genome.name + ".fasta")
    #
    #elif os.path.isfile(os.path.join(annotation_path, anno_genome.name + ".fa")):
    #    reference_location = os.path.join(annotation_path, anno_genome.name + ".fa")
    #
    #elif os.path.isfile(os.path.join(annotation_path, os.path.split(annotation_path)[1] + ".fasta")):
    #    reference_location = os.path.join(annotation_path, os.path.split(annotation_path)[1] + ".fasta")
    #
    #elif os.path.isfile(os.path.join(annotation_path, os.path.split(annotation_path)[1] + ".fa")):
    #    reference_location = os.path.join(annotation_path, os.path.split(annotation_path)[1] + ".fa")
    #
    #else:
    #    reference_location = os.path.join(annotation_path, os.path.split(annotation_path)[1] + ".fasta")
    #    print(reference_location)
    #    assert os.path.isfile(reference_location) # blast reference db files not found
    
    # print(ref)
    file_name = str(ref_position)
    query_len = int(query_len)

    query_seq = ref_genome[ref_position-1-query_len: ref_position-1]
    # print(query_seq)

    # MADE CHANGES HERE TO SAVE BLAST RESULTS AS TEMP FILES
    query_path = os.path.join(blast_dir_name, "temp_query.fasta")
    #query_path = f"{blast_dir_name}query_{file_name}.fasta"
    with open(query_path, "w") as o:
        SeqIO.write(query_seq, o, "fasta")
    output_path = os.path.join(blast_dir_name, "temp_results.xml")
    #output_path = f"{blast_dir_name}results_{file_name}.xml"
    
    os.system(
        f"blastn -db {reference_location} -query {query_path} -out {output_path} -outfmt 5")
    result_handle = open(output_path)
    blast_record = NCBIXML.read(result_handle)
    col_names = ['#HSPs', 'Percent_identity', 'Anno_position']
    n_alignments, identities, anno_position = verify_alignment(
        blast_record, query_len)

    # Does the base at the reference position match?
    if anno_position != 'Undetermined':
        anno_position = int(anno_position)
        if str(anno_genome[anno_position-1:anno_position-1+len(ref)].seq) != ref:
            return pd.Series({k: v for k, v in zip(col_names, ['Undetermined']*3)})
    # print(anno_position)
    return pd.Series({k: v for k, v in zip(col_names, [n_alignments, identities, anno_position])})


def verify_alignment(blast_record, query_len: int) -> Tuple[int, int, int]:
    try:
        n_hsps = 0
        assert len(blast_record.alignments) > 0
        # I think you can only have 0 or 1 alignments with a single reference. This is to catch any exceptions
        if len(blast_record.alignments) > 1:
            print('MORE THAN ONE ALIGNMENT!')
            raise ValueError
        n_hsps = len(blast_record.alignments[0].hsps)
        hsp = blast_record.alignments[0].hsps[0]
        assert hsp.align_length == query_len
        identities = hsp.identities*100/query_len
        return n_hsps, identities, hsp.sbjct_end + 1
    except AssertionError:
        return n_hsps, 'Undetermined', 'Undetermined'


def extract_features(row: pd.Series, selected_features_pr: pr.PyRanges) -> pd.Series:
    """Extracts relevant features (both within gene and intergenic) given a genomic position"""
    anno_position = row['Anno_position']
    if anno_position == 'Undetermined':
        return pd.Series({k: v for k, v in zip(['Location', 'Feature_type', 'Feature_name', 'Distance'], ['Undetermined']*4)})
    # selected_features and anno_position are 1-indexed
    anno_position = int(float(anno_position))
    pos_pr = pr.PyRanges(chromosomes=selected_features_pr.chromosomes[0],
                         starts=[anno_position],
                         ends=[anno_position])
    intersections = selected_features_pr.intersect(pos_pr).df
    if len(intersections) == 0:
        return extract_intergenic(selected_features_pr, pos_pr)
    intersections = intersections.drop(columns=['End', 'Chromosome'])
    if len(intersections) > 0:
        intersections['Location'] = 'Within_gene'
        intersections['Distance'] = np.nan
    if len(intersections) > 1:
        intersections = intersections.groupby(['Location', 'Start'])[
            'Feature_name', 'Feature_type', 'Distance'].agg(list).reset_index()
    return intersections[['Location', 'Feature_type', 'Feature_name', 'Distance']].T.squeeze()


def extract_intergenic(selected_features_pr: pr.PyRanges, pos_pr: pr.PyRanges) -> pd.Series:
    """Extracts information of neighbouring features given a genomic position"""
    #previous_gene = pos_pr.nearest(selected_features_pr, how='previous').df.drop(
    #    columns=['Chromosome', 'Start_b', 'End_b', 'End'])
    #next_gene = pos_pr.nearest(selected_features_pr, how='next').df.drop(
    #    columns=['Chromosome', 'Start_b', 'End_b', 'End'])
    previous_gene = pos_pr.nearest(selected_features_pr, how='previous').df
    next_gene = pos_pr.nearest(selected_features_pr, how='next').df
    if (len(previous_gene) > 0 ) and (len(next_gene) > 0 ):
        nearest = pd.concat([
            previous_gene[['Start', 'Feature_name', 'Feature_type', 'Strand', 'Distance']],
            next_gene[['Start', 'Feature_name', 'Feature_type', 'Strand', 'Distance']]], axis=0)

    elif len(previous_gene) > len(next_gene):
        nearest = previous_gene[['Start', 'Feature_name', 'Feature_type', 'Strand', 'Distance']]

    else:
        nearest = next_gene[['Start', 'Feature_name', 'Feature_type', 'Strand', 'Distance']]

    #print(nearest.columns)
    nearest['Location'] = 'Intergenic'
    nearest = nearest.groupby(['Location', 'Start'])[
        ['Feature_name', 'Feature_type', 'Distance']].agg(list).reset_index()
    #print(nearest)
    return nearest[['Location', 'Feature_type', 'Feature_name', 'Distance']].T.squeeze()


def annotate_aa(row: pd.Series, selected_features_pr: pr.PyRanges, anno_genome: SeqRecord) -> pd.Series:
    """For a given row in the output_df adds columns corresponding to details of the mutation"""
    aa_columns = {}
    if row['Anno_position'] == 'Undetermined':
        aa_columns['Mutation_type'] = 'Undetermined'
        aa_columns['AA_change'] = 'Undetermined'
        return pd.Series(aa_columns)
    if (row['Feature_type'] != 'CDS') or (row['Location'] == 'Intergenic'):
        aa_columns['Mutation_type'] = np.nan
        aa_columns['AA_change'] = np.nan
        return pd.Series(aa_columns)
    ref = row['Ref']
    alt = row['Alt']

    if ',' in alt:
        print("WARNING, unexpected character ',' detected at coordinate {}".format(row['Ref_pos']))
        print("Full Alt value: ", alt)
        print("Continuing forward with final entry in 'Alt' field")
        alt = alt.split(',')[-1]
    
    if len(ref) == len(alt):
        anno_position = int(float(row['Anno_position']))
        #print(row)
        ref_protein, mut_protein = mutate_protein(
            anno_position, ref, alt, row['Feature_name'], selected_features_pr, anno_genome)
        aa_columns['AA_change'], aa_columns['Mutation_type'] = extract_missense(
            ref_protein, mut_protein)
    if len(ref) < len(alt):
        aa_columns['AA_change'] = np.nan
        if len(ref)-len(alt) % 3 == 0:
            aa_columns['Mutation_type'] = 'In-frame Ins'
        else:
            aa_columns['Mutation_type'] = 'Frameshift Ins'
    if len(ref) > len(alt):
        aa_columns['AA_change'] = np.nan
        if len(ref)-len(alt) % 3 == 0:
            aa_columns['Mutation_type'] = 'In-frame Del'
        else:
            aa_columns['Mutation_type'] = 'Frameshift Del'
    return pd.Series(aa_columns)


def mutate_protein(anno_position: int, ref: str, alt: str, feature_name: str, selected_features_pr: pr.PyRanges, anno_genome: SeqRecord) -> Tuple[str, str]:
    """Using the annotated position, creates the WT and the mutant protein"""
    # first extract the feature
    selected_features = selected_features_pr.df
    feature_start = selected_features.loc[selected_features['Feature_name']
                                          == feature_name, 'Start'].values[0]
    feature_end = selected_features.loc[selected_features['Feature_name']
                                        == feature_name, 'End'].values[0]
    feature_strand = selected_features.loc[selected_features['Feature_name']
                                           == feature_name, 'Strand'].values[0]
    # the genomes are zero-indexed
    gene_seq = str(anno_genome[feature_start-1: feature_end-1].seq)
    anno_position_gene = anno_position-feature_start
    # Check if references match
    assert gene_seq[anno_position_gene:anno_position_gene+len(ref)] == ref
    mutated_gene_seq = gene_seq[:anno_position_gene] + \
        alt + gene_seq[anno_position_gene+len(ref):]
    if feature_strand == '-':
        gene_seq = str(Seq(gene_seq).reverse_complement())
        mutated_gene_seq = str(Seq(mutated_gene_seq).reverse_complement())
    ref_protein = SeqRecord(Seq(gene_seq)).translate().seq
    mut_protein = SeqRecord(Seq(mutated_gene_seq)).translate().seq
    return str(ref_protein), str(mut_protein)


def extract_missense(ref_protein: str, mut_protein: str) -> List:
    """For a WT and mutant protein identifies missense and silent mutations"""
    missense_aa = []
    for n in range(len(ref_protein)):
        ref_aa = ref_protein[n]
        mut_aa = mut_protein[n]
        if ref_aa != mut_aa:
            # make it 1 indexed
            missense_aa.append(f"{ref_aa}{n+1}{mut_aa}")
    if len(missense_aa) == 0:
        return [np.nan, 'Silent']
    return [','.join(missense_aa), 'Missense']


def add_lab_ref(row: pd.Series, lab_ref: pd.DataFrame) -> bool:
    """Adds information regarding the mutations presence in the lab_reference"""
    subset_row = row[['Ref_pos', 'Ref', 'Alt']]
    subset_lab_ref = lab_ref[['Ref_pos', 'Ref', 'Alt']]
    return subset_lab_ref.apply(lambda x: subset_row.equals(x), axis=1).any()


if __name__ == "__main__":
    main()

# def verify_alignment(blast_record, query_len: int) -> int:
#     """Verifies if alignment is a high-confidence one, and returns the position in the annotated genome"""
#     try:
#         assert len(blast_record.alignments) == 1
#         assert len(blast_record.alignments[0].hsps) == 1
#         hsp = blast_record.alignments[0].hsps[0]
#         assert hsp.identities == query_len
#         # print(hsp.query[0:100] + '...')
#         # print(hsp.match[0:100] + '...')
#         # print(hsp.sbjct[0:100] + '...')
#         # note this is 1-based
#         return hsp.sbjct_end + 1
#     except AssertionError:
#         return 'Undetermined'


# def annotate_mutations(row: pd.Series, ref_genome: SeqRecord, anno_genome: SeqRecord, selected_features_pr: pr.PyRanges, blast_dir_name: str) -> pd.Series:
#     """For a given mutation, calls BLAST against the annotated genome and extracts features of the identified position"""
#     ref = row['REF']
#     # Is ref base identical between annotated and unannotated?

#     extracted_features = extract_features(selected_features_pr, anno_position)
#     return extracted_features


# # run_vcf_to_csv('./vcf/accD6_diploid.lenient.vcf')
# ref_genome = SeqIO.read(
#     './../../../Reference/H37RvCO/H37RvCO.fasta', format="fasta")

# anno_genome = SeqIO.read(
#     './../../../Reference/AL123456/AL123456.gbk', "genbank")

# selected_features = pd.read_csv(
#     './../../../AZ-Scripts/metadata/AL123456_selected_features.csv')
# # print(selected_features['Feature_type'].value_counts())

# selected_features_pr = pr.PyRanges(selected_features)
# vcf = pd.read_csv('./vcf/accD6_haploid.stringent.vcf', sep='\t', header=28)
# vcf.head()

# vcf_len = pd.read_csv('./vcf/accD6_diploid.lenient.vcf', sep='\t', header=28)
# # a = vcf.apply(annotated_mutations, axis=1, args=(ref_genome, anno_genome, selected_features_pr))
