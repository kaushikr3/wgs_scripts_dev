from __future__ import annotations

import argparse
import os
import glob
from collections import defaultdict

import numpy as np
import pandas as pd
import pyranges as pr
from Bio import SearchIO, SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



def blast_h37rv(row: pd.Series, ref_genome: SeqRecord, anno_genome: SeqRecord, blast_dir_name: str, 
              reference_location: str, query_len=500) -> pd.Series:
    """For each row, uses the reference position to create a query sequence of length query length. 
    This is then blast-ed against the annotated genome. The alignment is verified and details are returned as a Series"""
    # Note that reference genome is zero-indexed and ref position is 1-indexed
    ref_position = row['Pos']
    ref = row['Ref']

    file_name = str(ref_position)
    query_len = int(query_len)
    half_query_len = query_len // 2

    query_seq = ref_genome[ref_position-1-half_query_len: ref_position-1+half_query_len]

    query_path = f"{blast_dir_name}query_{file_name}.fasta"
    with open(query_path, "w") as o:
        SeqIO.write(query_seq, o, "fasta")
    output_path = f"{blast_dir_name}results_{file_name}.xml"
    os.system(
        f"blastn -db {reference_location} -query {query_path} -out {output_path} -outfmt 5")
    result_handle = open(output_path)
    blast_record = NCBIXML.read(result_handle)
    col_names = ['#HSPs', 'Percent_identity', 'Anno_position']
    n_alignments, identities, anno_position = verify_homolog_alignment(
        blast_record, query_len)

    # Does the base at the reference position match?
#     if anno_position != 'Undetermined':
#         anno_position = int(anno_position)
#         if str(anno_genome[anno_position-1:anno_position-1+len(ref)].seq) != ref:
#             return pd.Series({k: v for k, v in zip(col_names, ['Undetermined']*3)})
    # print(anno_position)
    return pd.Series({k: v for k, v in zip(col_names, [n_alignments, identities, anno_position])})


def verify_homolog_alignment(blast_record, query_len: int) -> Tuple[int, int, int]:
    try:
        n_hsps = 0
        assert len(blast_record.alignments) > 0
        # I think you can only have 0 or 1 alignments with a single reference. This is to catch any exceptions
        if len(blast_record.alignments) > 1:
            print('MORE THAN ONE ALIGNMENT!')
            raise ValueError
        n_hsps = len(blast_record.alignments[0].hsps)
        hsp = blast_record.alignments[0].hsps[0]
        assert hsp.align_length > 0.9*query_len
        identities = hsp.identities*100/query_len
        return n_hsps, identities, hsp.sbjct_end + 1
    except AssertionError:
        return n_hsps, 'Undetermined', 'Undetermined'
    

def extract_features(row: pd.Series, selected_features_pr: pr.PyRanges, position_column: str) -> pd.Series:
    """Extracts relevant features (both within gene and intergenic) given a genomic position"""
    anno_position = row[position_column]
    
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
    previous_gene = pos_pr.nearest(selected_features_pr, how='previous').df.drop(
        columns=['Chromosome', 'Start_b', 'End_b', 'End'])
    next_gene = pos_pr.nearest(selected_features_pr, how='next').df.drop(
        columns=['Chromosome', 'Start_b', 'End_b', 'End'])
    nearest = pd.concat([previous_gene, next_gene], axis=0)
    nearest['Location'] = 'Intergenic'
    nearest = nearest.groupby(['Location', 'Start'])[
        ['Feature_name', 'Feature_type', 'Distance']].agg(list).reset_index()
    return nearest[['Location', 'Feature_type', 'Feature_name', 'Distance']].T.squeeze()


def add_h37rv_homolog_data(df: pd.DataFrame, 
                           ref_genome: SeqRecord, h37rv_genome: SeqRecord, 
                           RV_selected_features_pr: pr.PyRanges, 
                           reference_location: str, blast_dir_name: str) -> pd.DataFrame:   
    
    df[['H37Rv_homolog_HSPs', 'H37Rv_homolog_%identity', 'H37Rv_homolog_position']] = df.apply(
        blast_h37rv, axis=1, args=(ref_genome, h37rv_genome, blast_dir_name, reference_location))

    df[['H37Rv_Location', 'H37Rv_Feature_type', 'H37Rv_Feature_name', 'H37Rv_Distance']] = df.apply(
        extract_features, axis=1, args=(RV_selected_features_pr, 'H37Rv_homolog_position'))


    cols_to_explode = ['H37Rv_Feature_type', 'H37Rv_Feature_name', 'H37Rv_Distance']
    df = df.set_index([c for c in df.columns if c not in cols_to_explode]).apply(pd.Series.explode).reset_index()

    gene_details = pd.read_csv("~/wgs/metadata_wgs/AL123456_gene_details.csv")
    gene_details.columns = ["ORF", "H37Rv_Gene_Name", "H37Rv_Feature_Description"]
    df = pd.merge(df, gene_details, left_on=['H37Rv_Feature_name'], right_on=['ORF'], how='left').drop(columns='ORF')

    cols_to_group = ['H37Rv_Distance', 'H37Rv_Feature_type', 'H37Rv_Feature_name', 'H37Rv_Gene_Name', 'H37Rv_Feature_Description']

    df = df.fillna('NA')
    # group duplicates
    df = df.groupby([c for c in df.columns if c not in cols_to_group])[cols_to_group].agg(list).reset_index()
    df[cols_to_group] = df[cols_to_group].applymap(lambda x: ";".join(map(str, x)))
    
    # reorder columns:
    df = df[['Strain', 'Ref_genome', 'Ref_pos', 'Ref', 'Alt', '#Ref_reads','#Alt_reads', 'Coverage', 'Alt_freq', 'GQ', 
             'Anno_genome', '#HSPs', 'Percent_identity', 'Anno_position', 'Location', 'Distance',
             'Feature_type', 'Feature_name', 'Feature_gene_name', 'Feature_description', 'AA_change', 'Mutation_type', 
             'H37Rv_homolog_HSPs', 'H37Rv_homolog_%identity', 'H37Rv_homolog_position', 'H37Rv_Location', 
             'H37Rv_Distance', 'H37Rv_Feature_type', 'H37Rv_Feature_name','H37Rv_Gene_Name', 'H37Rv_Feature_Description',
             'Lab_reference','Present_in_lab_reference_stringent','Present_in_lab_reference_lenient']]
    
    return df


def main():
    parse = argparse.ArgumentParser(prog='H37Rv_homologs')
    
    parse.add_argument('-ref_fa', required=True, dest='ref_path', type=str,
                       help='Path to reference fasta')
    
    parse.add_argument('-mtb_gbk', required=True, dest='h37rv_gb', type=str,
                       help='Path to H37Rv Genbank file')
    
    parse.add_argument('-csv', required=True, dest='csv_dir', type=str,
                       help='Path to current csv files')
    
    parse.add_argument('-csv_new', required=True, dest='csv_dir_new', type=str,
                       help='Path to directory to write new csv files into')
    
    parse.add_argument('-blast_new', required=True, dest='blast_dir_new', type=str,
                       help='Path to directory to write H37Rv blast files into')
    
    args = parse.parse_args()
    
    # read in annotation data:
    ref_genome = SeqIO.read(args.ref_path, format="fasta")
    reference_path = os.path.split(os.path.split(args.ref_path)[0])[0]

    h37rv_genome = SeqIO.read(args.h37rv_gb, format="genbank")
    reference_location = os.path.join(os.path.split(args.h37rv_gb)[0], h37rv_genome.name + ".fasta")


    RV_selected_features = pd.read_csv("~/wgs/metadata_wgs/AL123456_selected_features.csv")
    h37rv_feature_pr = pr.PyRanges(RV_selected_features)
    
    # make new csv and blast directories
    os.system(f"mkdir {args.csv_dir_new}")
    os.system(f"mkdir {args.blast_dir_new}")

    csvs = [f for f in os.listdir(args.csv_dir) if os.path.isfile(os.path.join(args.csv_dir, f))]

    for file in csvs:
        blast_dir = os.path.join(args.blast_dir_new, "_".join(file.split('.')[:-1]))
        os.system("mkdir {}".format(blast_dir))
        
        print("Running {}".format(file))
        df = pd.read_csv(os.path.join(args.csv_dir, file))
        df = add_h37rv_homolog_data(df, 
                                    ref_genome, h37rv_genome, 
                                    h37rv_feature_pr, 
                                    reference_location, blast_dir)

        df.to_csv(os.path.join(args.csv_dir_new, file), index=False)


if __name__ == "__main__":
    main()
