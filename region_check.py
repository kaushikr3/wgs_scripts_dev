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


# note that the plasmid insertion files have locations that are 1-indexed with the end as inclusive

def main():
    my_parser = argparse.ArgumentParser(prog='region_check')
    my_parser.add_argument('-regions_file', required=True,
                           dest='regions_file', type=str)
    my_parser.add_argument('-samples_file', required=True,
                           dest='samples_file', type=str)
    my_parser.add_argument('-stringency', required=True,
                           dest='stringency', choices={'lenient', 'stringent'}, type=str)
    my_parser.add_argument('-cov_dir', required=False,
                           default='cov', dest='cov_dir', type=str)
    my_parser.add_argument('-csv_dir', required=False,
                           default='csv', dest='csv_dir', type=str)
    args = my_parser.parse_args()

    regions = pd.read_csv(args.regions_file)
    samples = pd.read_csv(args.samples_file)
    out_df = pd.merge(samples, regions, how='left', on='Strain_name')
    out_df['#Low_coverage'] = out_df.apply(
        check_coverage, axis=1, args=(args.cov_dir,))
    snp_dfs = [check_snps(row, args.csv_dir, args.stringency) for _, row in out_df.iterrows()]
    all_snps_df = pd.concat(snp_dfs)
    # print(all_snps_df)
    num_snps = all_snps_df.groupby(['Strain', 'Inserted_plasmid']).size().rename('#Mutations').reset_index()
    # print(num_snps)
    num_snps = num_snps.rename(columns={'Strain': 'Sample_name'})
    out_df = pd.merge(out_df, num_snps, how='left', on=['Sample_name', 'Inserted_plasmid'])
    result_to_excel(f'Plasmid_status_{args.stringency}.xlsx', [
                    out_df, all_snps_df], ['Summary', 'Details'])
    return None


def check_coverage(row, cov_dir):
    """For a given sample, reads the corresponding coverage file, and returns the number of sites with <5 coverage in a given range"""
    sample_name = row['Sample_name']
    cov = pd.read_csv(f'./{cov_dir}/{sample_name}.cov', header=None,
                      names=['Strain', 'Location', 'Coverage'], sep='\t')
    start = row['Start_minus100']
    end = row['End_plus100']
    cov_filt = cov[(cov['Location'] >= start) & (cov['Location'] <= end)]
    return (cov_filt['Coverage'] < 5).sum()


def check_snps(row, csv_dir, stringency):
    """For a given sample, reads the corresponding csv file and returns the filtered snp df"""
    sample_name = row['Sample_name']
    if stringency == 'stringent':
        csv = pd.read_csv(
            f'./{csv_dir}/{sample_name}_haploid.stringent.annotated.csv')
    elif stringency == 'lenient':
        csv = pd.read_csv(
            f'./{csv_dir}/{sample_name}_diploid.lenient.annotated.csv')
    start = row['Start_minus100']
    end = row['End_plus100']
    snp_filt = csv[(csv['Ref_pos'] >= start) & (csv['Ref_pos'] <= end)]
    snp_filt['Inserted_plasmid'] = row['Inserted_plasmid']
    snp_filt['Site'] = row['Site']
    snp_filt['Plasmid_start'] = row['Start']
    snp_filt['Plasmid_end'] = row['End']
    snp_filt['Position_in_plasmid'] = snp_filt['Ref_pos'] - snp_filt['Plasmid_start'] +1 
    snp_filt = snp_filt[['Strain', 'Ref_genome', 'Inserted_plasmid', 'Site', 'Plasmid_start', 'Plasmid_end', 'Ref_pos', 'Position_in_plasmid', 'Ref', 'Alt', '#Ref_reads', '#Alt_reads', 'Alt_freq', 'GQ']]
    return snp_filt


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
    writerReport.close()
    return output_name


if __name__ == "__main__":
    main()
