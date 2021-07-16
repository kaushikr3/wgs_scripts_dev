import os
import argparse
import pandas as pd
import numpy as np

# VCF Column formatting:
format_cols_dict = {
    "CHROM": "Chrom",
    "POS": "Pos",
    "REF": "Ref",
    "ALT": "Alt",
    "QUAL": "Qual",
    "sample.GT": "GT",
    "sample.GQ": "GQ",
    "sample.PL": "PL",
    "sample.AD": "AD",
    "sample.DP": "DP",
    "snippy.GT": "GT",
    "snippy.PL": "PL",
    "snippy.DP": "DP",  # depth of coverage
    "snippy.AO": "AC",  # alternate allele read count
    "snippy.RO": "RC"  # reference allele read count
}

### NEED TO PULL IN COORDINATES OF INTEREST; THATS ALL

def main():
    """

    :return: Writes excel sheet holding annotated SNP data from parsed VCF files
    """
    parse = argparse.ArgumentParser(prog='Recombinant_construct_SNP_annotation')

    parse.add_argument('-start', required=True, dest='start', type=str,
                       help='Start coordinate of recombinant construct, or list of start coordinates for multiple'
                            'constructs separated by commas -- list order must correspond to order of end coords')

    parse.add_argument('-end', required=True, dest='end', type=str,
                       help='End coordinate of recombinant construct, or list of end coordinates for multiple'
                            'constructs separated by commas -- list order must correspond to order of start coords')

    parse.add_argument('-vcf', required=True, dest='vcf_dir', type=str,
                       help='Path to directory holding PARSED vcfs')

    parse.add_argument('-out', required=True, dest='out', type=str,
                       help='Path to directory to write excel files into')

    args = parse.parse_args()
    os.system(f"mkdir {args.out}")

    start = [int(coord) for coord in args.start.split(',')]
    end = [int(coord) for coord in args.end.split(',')]

    assert len(start) == len(end) # Must provide same number of start and end coordinates

    files = [f for f in os.listdir(args.vcf_dir) if os.path.isfile(os.path.join(args.vcf_dir, f))]
    samples = pd.Series(files).str.split('.', expand=True
                                         )[0].str.rstrip('freebayes').str.rstrip('gatk').str.rstrip('_').unique()

    for sample in samples:
        print("Running sample: ", sample)
        vcf_dict = {}

        if args.h37rv_homology or args.ref_strain == 'H37Rv':
            blast_sub_dir = os.path.join(args.blast_dir, sample)
            os.system(f"mkdir {blast_sub_dir}")

        for file in files:
            if sample in file:
                if 'haploid' in file:
                    vcf_dict.update({"gatk_haploid": os.path.join(args.vcf_dir, file)})
                elif 'diploid' in file:
                    vcf_dict.update({"gatk_diploid": os.path.join(args.vcf_dir, file)})
                elif 'freebayes' in file:
                    vcf_dict.update({"freebayes": os.path.join(args.vcf_dir, file)})

        print(vcf_dict)
        df = get_full_annotated_csv(vcf_dict['gatk_haploid'], vcf_dict['gatk_diploid'], vcf_dict['freebayes'],
                                    start, end)

        write_excel_file(df, os.path.join(args.out, sample + "_SNP.xlsx"))


def get_full_annotated_csv(gatk_haploid_parsed_path, gatk_diploid_parsed_path, freebayes_parsed_path,
                           start, end):
    """
    Generates fully annotated SNP df with reference annotation, H37Rv homology, and lab strain intersection

    :param gatk_haploid_parsed_path: Path to gatk haploid/stringent parsed vcf file
    :param gatk_diploid_parsed_path: Path to gatk diplid/lenient parsed vcf file
    :param freebayes_parsed_path: Path to freebayes parsed vcf file
    :param start: list of coordinates where recombinant regions begin [int, ....., int]
    :param end: list of coordinates where recombinant regions end [int, ....., int]
    :return: Fully annotated SNP df, ready to be written to outfile
    """

    # open SNPs
    hgatk = pd.read_csv(gatk_haploid_parsed_path, delimiter='\t')
    dgatk = pd.read_csv(gatk_diploid_parsed_path, delimiter='\t')
    free = pd.read_csv(freebayes_parsed_path, delimiter='\t')

    # add SNP_tool, accuracy and error columns, and rename columns
    hgatk = format_vcf_fields(hgatk, 'gatk')
    dgatk = format_vcf_fields(dgatk, 'gatk')
    free = format_vcf_fields(free, 'freebayes')

    # merge SNP sets from different sources together:
    df = merge_vcf_dfs(hgatk, dgatk, free)

    # reorder columns and drop composite_hits:
    df = df[['Chrom', 'Pos', 'Ref', 'Alt', 'Qual', 'Error_prob', 'Accuracy_prob', 'Tool',
             'GT', 'GQ', 'PL', 'DP', 'RC', 'AC']]

    # limit SNP rows to those inside the recombinant construct coordinates
    conds = []
    for i in range(len(start)):
        conds.append((df['Pos'] > start[i]) & (df['Pos'] < end[i]))

    df = pd.concat([df[conds[i]] for i in range(len(start))])

    df = df.fillna('NA')

    return df


def merge_vcf_dfs(haploid_gatk, diploid_gatk, free):
    """

    :param haploid_gatk: Parsed VCF df from haploid GATK SNP caller
    :param diploid_gatk: Parsed VCF df from haploid GATK SNP caller
    :param free: Parsed VCF df from Freebayes SNP caller
    :return: Df with SNPs from varied callers merged, and duplicated SNPs labeled, with duplicate entries dropped
    """

    # merge GATK SNP dfs first, and drop lenient version of duplicates, then merge in freebayes SNPs
    df = haploid_gatk.append(diploid_gatk).reset_index(drop=True)
    df = df.sort_values(by=['Tool', 'Pos'])[~df[['Chrom', 'Pos', 'Ref', 'Alt']].duplicated()]
    df = df.append(free).reset_index(drop=True)

    # label all duplicate SNPs and drop Freebayes version of duplicates
    df['Composite_hit'] = df[['Chrom', 'Pos', 'Ref', 'Alt']].duplicated(keep=False)
    df = df[~(df['Composite_hit'] & (df['Tool'] == 'Freebayes'))]
    # df = df[~((df['Composite_hit'] == True) & (df['Tool'] == 'Freebayes'))]

    # Alter tool column to reflect duplicate SNP detections
    df['Tool'] = df.apply(lambda x:
                          x['Tool'] + '+Freebayes' if x['Composite_hit'] is True else x['Tool'], axis=1)

    return df


def get_error_probability(quality):
    """
    :param quality: phred-scaled estimate of confidence that a polymorphism exists at a given genome position
    :return: probability that the (variant) base has been called incorrectly by the sequencer
    """
    error = 10 ** (-quality / 10)

    return error


def format_vcf_fields(df, source):
    if source == 'gatk':
        column_check = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'sample.GT', 'sample.GQ', 'sample.PL', 'sample.AD',
                        'sample.DP']
        assert df.columns.tolist() == column_check  # vcf df in unexpected format

        # split allele depth column into reference allele count and alternate allele count columns
        df[["RC", "AC"]] = df['sample.AD'].str.split(',', expand=True)[[0, 1]]
        df = df.drop('sample.AD', axis=1)

        df['Tool'] = 'GATK'

    elif source == 'freebayes':
        column_check = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'sample.GT', 'sample.GQ', 'sample.PL', 'sample.DP',
                        'sample.AD']
        assert df.columns.tolist() == column_check  # vcf df in unexpected format

        # split allele depth column into reference allele count and alternate allele count columns
        df[["RC", "AC"]] = df['sample.AD'].str.split(',', expand=True)[[0, 1]]
        df = df.drop('sample.AD', axis=1)

        df['Tool'] = 'Freebayes'

    elif source == 'snippy':
        column_check = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'snippy.GT', 'snippy.PL', 'snippy.DP', 'snippy.RO',
                        'snippy.AO']
        assert df.columns.tolist() == column_check  # vcf df in unexpected format

        df['Tool'] = 'Snippy'

    df = df.rename(columns=format_cols_dict)

    df['error_prob'] = df['Qual'].apply(lambda x: get_error_probability(x))
    df['accuracy_prob'] = 1 - df['error_prob']

    # returns vcf tsvs with common sense column names, and error and accuracy metric cols added in
    # any other columns worth adding??

    return df


def write_excel_file(df, filename):
    """
    Writes SNPs to Excel file with 5 sheets holding 5 levels of confidence that SNP was called correctly
    :param df: fully annotated SNP df ready to be written to Excel file (pd.DataFrame)
    :param filename: filename dest (str)
    :return: None
    """

    dfs_to_write = {
        'Accuracy Prob > 0.75': df[df['accuracy_prob'] > 0.75].sort_values(by='Pos'),
        'Accuracy Prob > 0.50': df[df['accuracy_prob'] > 0.50].sort_values(by='Pos'),
        'Accuracy Prob > 0.25': df[df['accuracy_prob'] > 0.25].sort_values(by='Pos'),
        'Accuracy Prob > 0.05': df[df['accuracy_prob'] > 0.05].sort_values(by='Pos'),
        'All Possible SNPs': df.sort_values(by='Pos'),
                   }

    writer = pd.ExcelWriter(filename)
    # loop through drug's dfs
    for sheet, df in dfs_to_write.items():
        df.to_excel(writer, sheet_name=sheet, index=False)

        # set up df writer settings
        worksheet = writer.sheets[sheet]  # pull worksheet object
        for idx, col in enumerate(df):
            worksheet.set_column(idx, idx, 14)  # set column width

    writer.save()
    print("Excel File Written")


if __name__ == "__main__":
    main()

