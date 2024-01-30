from __future__ import annotations
import pandas as pd
import argparse
from pathlib import Path
import glob
import pyranges as pr
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from Bio import SeqIO
pd.set_option('mode.chained_assignment', None)


# note that feature range files but not pdim_location files have to be +1 at the end.

def main():
    my_parser = argparse.ArgumentParser(prog='pdim_check')
    my_parser.add_argument('-strain', required=True, dest='strain', 
                           choices={'H37RvCO', 'HN878', 'Erdman', 'BCG', 'AL123456', 'NC_000962'},
                           type=str)
    my_parser.add_argument('-stringency', required=True,
                           dest='stringency', choices={'Lenient', 'Stringent'}, type=str)
    my_parser.add_argument('-cov_dir', required=False,
                           default='cov', dest='cov_dir', type=str)
    my_parser.add_argument('-csv_dir', required=False,
                           default='csv', dest='csv_dir', type=str)
    my_parser.add_argument('-genome', required=True,
                           dest='genome_path', type=str)

    args = my_parser.parse_args()

    names = [Path(x).stem for x in glob.glob(f'./{args.cov_dir}/*.cov')]
    pdim_range = pd.read_csv('~/metadata_wgs/PDIM_locations_new.csv')
    covs = [pd.read_csv(f'{args.cov_dir}/{name}.cov', header=None,
                        names=['Strain', 'Location', 'Coverage'], sep='\t') for name in names]

    csvs = [pd.read_excel(f'{args.csv_dir}/{name}_SNP.xlsx', sheet_name=args.stringency) for name in names]
    csvs_dict = dict(zip(names, csvs))

    selected_features = pd.read_csv(
        f"~/metadata_wgs/PDIM_{args.strain}_features.csv")
    selected_features_pr = pr.PyRanges(selected_features)
    genome = SeqIO.read(args.genome_path, format="genbank")

    start = pdim_range.loc[pdim_range['Name']
                           == args.strain, 'Start'].values[0]
    end = pdim_range.loc[pdim_range['Name'] == args.strain, 'End'].values[0]

    # for summary
    coverage = [check_coverage(cov, start, end) for cov in covs]
    snps = [check_snps(csv, start, end) for csv in csvs]

    print(names)
    print(coverage)
    print(snps)

    summary_df = pd.DataFrame(
        {'Strain': names, '#Low_coverage': coverage, '#Mutations': snps})
    print(type(summary_df))
    print(summary_df.shape)

    summary_df['PDIM_status'] = summary_df.apply(check_pdim, axis=1)

    # for snps
    snps_df = pd.concat([extract_snp_details(csv, csvs_dict[csv], start, end) for csv in csvs_dict.keys()])
    #snps_df = pd.concat([extract_snp_details(csv, start, end) for csv in csvs])
    if snps_df.shape[0] == 0:
        result_to_excel(f'PDIM_status_{args.stringency}.xlsx', [summary_df], [
                    'Summary'])
        return None
    snps_df[['Location', 'Feature_type', 'Feature_name', 'Distance']] = snps_df.apply(
        extract_features, axis=1, args=(selected_features_pr,))
    cols_to_explode = ['Feature_type', 'Feature_name', 'Distance']
    snps_df = snps_df.set_index([c for c in snps_df.columns if c not in cols_to_explode]).apply(
        pd.Series.explode).reset_index()
    snps_df[['AA_change', 'Mutation_type']
            ] = snps_df.apply(annotate_aa, axis=1, args=(selected_features_pr, genome))
    snps_df['Distance'] = snps_df['Distance'].astype(str)
    # add metadata
    snps_df = snps_df.fillna('NA')
    # group duplicates
    cols_to_group = ['Distance', 'Feature_type',
                     'Feature_name', 'AA_change', 'Mutation_type']
    snps_df = snps_df.groupby([c for c in snps_df.columns if c not in cols_to_group])[
        cols_to_group].agg(list).reset_index()
    snps_df[cols_to_group] = snps_df[cols_to_group].applymap(';'.join)
    cols_to_move = ['Present_in_lab_reference_stringent',
                    'Present_in_lab_reference_lenient']
    snps_df = snps_df[[
        col for col in snps_df.columns if col not in cols_to_move]+cols_to_move]

    result_to_excel(f'PDIM_status_{args.stringency}.xlsx', [summary_df, snps_df], [
                    'Summary', 'Details'])
    # writer = pd.ExcelWriter('PDIM_status.xlsx', engine='xlsxwriter')
    # summary_df.to_excel(writer, sheet_name='Summary', index=False)
    # snps_df.to_excel(writer, sheet_name='Details', index=False)
    # writer.save()


def check_coverage(cov, start, end):
    cov_pdim = cov[(cov['Location'] >= start) & (cov['Location'] <= end)]
    #return (cov_pdim['Coverage'] < 0).sum()
    return (cov_pdim['Coverage'] < 5).sum()


def check_snps(csv, start, end):
    snp_pdim = csv[(csv['Pos'] >= start) & (csv['Pos'] <= end)]
    return snp_pdim.shape[0]


def check_pdim(row):
    if (row['#Low_coverage'] == 0) and (row['#Mutations'] == 0):
        return 'Intact'
    else:
        return 'Dubious'


#def extract_snp_details(csv, start, end):
def extract_snp_details(sample_name, csv, start, end):
    csv['Strain'] = sample_name
    snp_pdim = csv[(csv['Pos'] >= start) & (csv['Pos'] <= end)]
    snp_pdim['#Mutations'] = snp_pdim.shape[0]
    lab_ref_cols = ['Present_in_lab_reference_stringent',
                    'Present_in_lab_reference_lenient']
    for col in lab_ref_cols:
        if col not in snp_pdim.columns:
            snp_pdim[col] = np.nan

    snp_pdim[lab_ref_cols] = snp_pdim[lab_ref_cols].astype(str)
    snp_pdim = snp_pdim[['Strain', '#Mutations', 'Ref_genome', 'Pos',
                         'Ref', 'Alt', 'Ref_reads', 'Alt_reads', 'GQ'] + lab_ref_cols]
    return snp_pdim


def extract_features(row: pd.Series, selected_features_pr: pr.PyRanges) -> pd.Series:
    """Extracts relevant features (both within gene and intergenic) given a genomic position"""
    ref_position = row['Pos']
    # selected_features and anno_position are 1-indexed
    pos_pr = pr.PyRanges(chromosomes=selected_features_pr.chromosomes[0],
                         starts=[ref_position],
                         ends=[ref_position])
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
    previous_gene = pos_pr.nearest(selected_features_pr, how='previous').df
    next_gene = pos_pr.nearest(selected_features_pr, how='next').df 

    nearest = pd.concat([previous_gene, next_gene], axis=0).drop(
        columns=['Chromosome', 'Start_b', 'End_b', 'End'])
    nearest['Location'] = 'Intergenic'
    nearest = nearest.groupby(['Location', 'Start'])[
        ['Feature_name', 'Feature_type', 'Distance']].agg(list).reset_index()
    return nearest[['Location', 'Feature_type', 'Feature_name', 'Distance']].T.squeeze()


def annotate_aa(row: pd.Series, selected_features_pr: pr.PyRanges, genome: SeqRecord) -> pd.Series:
    """For a given row in the snps_df adds columns corresponding to details of the mutation"""
    aa_columns = {}
    if (row['Feature_type'] != 'CDS') or ('Intergenic' in row['Feature_type']):
        aa_columns['Mutation_type'] = np.nan
        aa_columns['AA_change'] = np.nan
        return pd.Series(aa_columns)
    ref = row['Ref']
    alt = row['Alt']
    if len(ref) == len(alt):
        ref_position = int(row['Pos'])
        ref_protein, mut_protein = mutate_protein(
            ref_position, ref, alt, row['Feature_name'], selected_features_pr, genome)
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


def mutate_protein(ref_position: int, ref: str, alt: str, feature_name: str, selected_features_pr: pr.PyRanges, genome: SeqRecord) -> Tuple[str, str]:
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
    gene_seq = str(genome[feature_start-1: feature_end-1].seq)
    ref_position_gene = ref_position-feature_start
    mutated_gene_seq = gene_seq[:ref_position_gene] + \
        alt + gene_seq[ref_position_gene+len(ref):]
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


def result_to_excel(output_name, dataframes_list, sheet_names_list):
    out_path = output_name
    writerReport = pd.ExcelWriter(out_path, engine='xlsxwriter',
                                  datetime_format='yyyymmdd', date_format='yyyymmdd')
    workbook = writerReport.book
    # loop through the list of dataframes to save every dataframe into a new sheet in the excel file
    for i, dataframe in enumerate(dataframes_list):
        # choose the sheet name from sheet_names_list
        sheet_name = sheet_names_list[i]
        dataframe.to_excel(writerReport, sheet_name=sheet_name,
                           index=False, startrow=0)
        # Add a header format.
        format = workbook.add_format({'bold': True})
        # Write the column headers with the defined format.
        worksheet = writerReport.sheets[sheet_name]
        for col_num, col_name in enumerate(dataframe.columns.values):
            worksheet.write(0, col_num, col_name, format)
        worksheet.freeze_panes(1, 0)
        # loop through the columns in the dataframe to get the width of the column
        for j, col in enumerate(dataframe.columns):
            max_width = max([len(str(s))
                             for s in dataframe[col].values] + [len(col) + 2])
            # define a max width to not get to wide column
            if max_width > 50:
                max_width = 50
            worksheet.set_column(j, j, max_width)
    writerReport.save()
    return output_name


if __name__ == "__main__":
    main()
