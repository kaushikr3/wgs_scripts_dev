import os
import pandas as pd
import numpy as np
import pyranges as pr

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import GenBank

from generate_h37rv_homo



# paths to references files:
genbank_dict = {'H37Rv': 'reference_files/AL123456.gbk',
                # 'Erdman':
                'Msmeg': 'reference_files/NC_008596.gb',
                'BCG': 'reference_files/BCG_Pasteur.gb'}

fasta_dict = {'H37Rv': 'reference_files/H37RvCO.fasta',
              # 'Erdman':
              'Msmeg': 'reference_files/NC_008596.fa',
              'BCG': 'reference_files/BCG_Pasteur.fa'}

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


def main(gatk_parsed_path, freebayes_parsed_path, reference_strain, lab_variant_csv_path):
    gatk = pd.read_csv(gatk_parsed_path, delimiter='\t')
    free = pd.read_csv(freebayes_parsed_path, delimiter='\t')

    # add SNP_tool, accuracy and error columns, rename columns, and merge
    gatk = format_vcf_fields(gatk, 'gatk')
    free = format_vcf_fields(free, 'freebayes')
    df = gatk.append(free).reset_index(drop=True)

    # label all duplicate SNPs and drop Freebayes version of duplicates
    df['Composite_hit'] = df[['Chrom', 'Pos', 'Ref', 'Alt']].duplicated(keep=False)
    df = df[~((df['Composite_hit'] == True) & (df['Tool'] == 'Freebayes'))]
    df.loc[df.Composite_hit == True, "Tool"] = "GATK+Freebayes"

    df = df.drop('Composite_hit', axis=1)

    # pull in annotation data:
    genome = SeqIO.read(genbank_dict[reference_strain], "genbank")
    genome_seq = SeqIO.read(fasta_dict[reference_strain], format="fasta")

    annotation = get_annotation_df(genome)
    anno_pr = pr.PyRanges(annotation)
    translated_pr = pr.PyRanges(annotation[
                                    annotation['Feature_type'] != 'regulatory'][
                                    annotation['Feature_type'] != 'misc_feature'])

    # if reference is H37Rv, get H37RvCO to H37Rv homolog information and merge in annotation data:
    if reference_strain == 'H37Rv':


    # otherwise, annotate with given reference and generate in H37Rv homology data
    else:
        df = snp_df_annotation(df, genome, anno_pr, translated_pr)



    # # add lab_reference_info
    if lab_variant_csv_path is not None:
        lab_stringent_path = "reference_files/BCG_Lab_ref_csvs/20_S20_L001_haploid.stringent.annotated.csv"
        lab_lenient_path =  "reference_files/BCG_Lab_ref_csvs/20_S20_L001_diploid.lenient.annotated.csv"

        df = add_lab_ref(df, lab_stringent_path, lab_lenient_path)


    # sort SNPs by confidence and write into file with sheets thresholded by confidence


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


def get_annotation_df(genome: SeqRecord) -> pd.DataFrame:
    """Extracts relevant features from the genome with location and annotation data"""
    selected_features = []

    for feature in genome.features:
        feature_dict = {'Chromosome': genome.name}

        if feature.type == 'CDS':
            feature_dict['Feature_tag'] = feature.qualifiers['locus_tag'][0]

            # add some extra annotation -- standard BCG_gene#
            if 'old_locus_tag' in feature.qualifiers:
                feature_dict['Old_tag'] = feature.qualifiers['old_locus_tag'][0]

            # add some extra annotation -- gene name
            if 'gene' in feature.qualifiers:
                feature_dict['Gene_name'] = feature.qualifiers['gene'][0]

            # add some extra annotation -- gene description
            if 'product' in feature.qualifiers:
                feature_dict['Feature_description'] = feature.qualifiers['product'][0]

        elif feature.type == 'tRNA':
            if 'locus_tag' in feature.qualifiers:
                feature_dict["Feature_tag"] = feature.qualifiers['locus_tag'][0]

            else:
                feature_dict['Feature_tag'] = feature.qualifiers['product'][0]

            feature_dict['Gene_name'] = feature.qualifiers['product'][0]

            if 'product' in feature.qualifiers:
                feature_dict['Feature_description'] = feature.qualifiers['product'][0]

        elif (feature.type == 'ncRNA') or (feature.type == 'misc_RNA'):

            # get feature name
            if 'locus_tag' in feature.qualifiers:
                feature_dict['Feature_tag'] = feature.qualifiers['locus_tag'][0]

            else:
                feature_dict['Feature_tag'] = feature.qualifiers['gene'][0]

            # add some extra annotation -- feature name
            if 'gene' in feature.qualifiers:
                feature_dict['Gene_name'] = feature.qualifiers['gene'][0]

                # add some extra annotation -- feature description
            if 'product' in feature.qualifiers:
                feature_dict['Feature_description'] = feature.qualifiers['product'][0]

        elif feature.type == 'misc_feature':
            feature_dict['Feature_tag'] = feature.qualifiers['locus_tag'][0]

            if 'old_locus_tag' in feature.qualifiers:
                feature_dict['Old_tag'] = feature.qualifiers['old_locus_tag'][0]

            if 'note' in feature.qualifiers:
                feature_dict['Feature_description'] = feature.qualifiers['note'][0].split(";")[0]

        elif feature.type == 'regulatory':
            feature_dict['Feature_tag'] = feature.qualifiers['regulatory_class'][0]

            if 'note' in feature.qualifiers:
                feature_dict['Feature_description'] = feature.qualifiers['note'][0].split(";")[0]

        else:
            continue

        feature_dict['Feature_type'] = feature.type
        feature_dict['Start'] = int(feature.location.start)
        feature_dict['End'] = int(feature.location.end)

        if feature.strand == 1:
            feature_dict['Strand'] = '+'
        elif feature.strand == -1:
            feature_dict['Strand'] = '-'

        selected_features.append(feature_dict)

    feats = pd.DataFrame(selected_features)

    # merge old/new locus tags into feature_name column:
    if 'Old_tag' in feats:
        feats['Feature_name'] = feats.apply(lambda x: x['Feature_tag'] if pd.isna(x['Old_tag']) else x['Old_tag'],
                                            axis=1)

    # make this 1-indexed
    feats['Start'] = feats['Start'] + 1
    feats['End'] = feats['End'] + 1

    return feats


def snp_df_annotation(df, genbank, anno_pr, translated_pr):
    # pull in snp annotation data:
    df[['Feature_type', 'Distance_to_feature', 'Feature_name', 'Gene_name', 'Feature_description']
    ] = df.apply(lambda x: annotate_snp(x, anno_pr, translated_pr), axis=1)

    # seperate out data so that any snp with multiple features can be properly annotated
    cols_to_explode = ['Distance_to_feature', 'Feature_name', 'Gene_name', 'Feature_description']
    df = pd.DataFrame(df.set_index([c for c in df.columns if c not in cols_to_explode])
                      ).apply(pd.Series.explode).reset_index()

    # annotate mutation types in coding sequences
    df[['Mutation_type', 'AA_change']] = df.apply(lambda x: annotate_aa(x, anno_pr, genbank), axis=1)

    # merge annotation data back into a single line / snp
    cols_to_group = ['Feature_type', 'Distance_to_feature', 'Feature_name', 'Gene_name', 'Feature_description',
                     'Mutation_type', 'AA_change']
    df['Distance_to_feature'] = df['Distance_to_feature'].astype(str)
    df = df.fillna('NA')

    df = df.groupby([col for col in df.columns if col not in cols_to_group]).agg(list).reset_index()

    # join list entries in df
    df[cols_to_group] = df[cols_to_group].apply(lambda x: x.str.join(";"))

    # simplify feature_type column; drop {};{} field unless there are actually two different feature types
    df['Feature_type'] = df['Feature_type'].str.split(';', expand=True).apply(
        lambda x: x[0] if (x[0] == x[1] or pd.isna(x[1]))
        else "{};{}".format(x[0], x[1]), axis=1)

    return df


def annotate_snp(snp, anno_pr, translated_pr):
    snp_coord = snp['Pos']
    snp_pr = pr.PyRanges(chromosomes=anno_pr.chromosomes[0],
                         starts=[snp_coord], ends=[snp_coord])
    intersect = anno_pr.intersect(snp_pr).df

    # return flanking features if SNP is intergenic or occurs in a regulatory feature
    if len(intersect) == 0:
        snp_annotation = get_adjacent_features(anno_pr, snp_pr)
        snp_annotation = pd.concat([pd.Series({'Feature_type': ['Intergenic', 'Intergenic']}),
                                    snp_annotation])

        return snp_annotation

    # if only intersecting features are regulatory or misc_features, return nearest features instead
    elif set(intersect['Feature_type']).union({'regulatory', 'misc_feature'}) == {'misc_feature', 'regulatory'}:
        snp_annotation = get_adjacent_features(translated_pr, snp_pr)
        snp_annotation = pd.concat([pd.Series({'Feature_type':
                                                   ["Intergenic ({})".format(i) for i in
                                                    intersect['Feature_description']]}),
                                    snp_annotation])

        return snp_annotation

    # if one overlapping feature return that feature's details
    elif len(intersect) == 1:
        intersect['Distance'] = np.nan
        snp_annotation = intersect[
            ['Feature_type', 'Distance', 'Feature_name', 'Gene_name', 'Feature_description']].T.squeeze()

        return snp_annotation

    # else more than one feature; squeeze into one row
    else:
        intersect['Distance'] = np.nan
        snp_annotation = intersect.groupby(['Start'])[
            ['Feature_type', 'Distance', 'Feature_name', 'Gene_name', 'Feature_description']
        ].agg(list).reset_index().T.squeeze().drop('Start')

        return snp_annotation


def get_adjacent_features(translated_pr: pr.PyRanges, snp_pr: pr.PyRanges) -> pd.Series:
    """Extracts information of neighbouring features given an intergenic genomic position"""

    # pull annotation data on previous and next features:
    adjacent_feats = pd.concat([snp_pr.nearest(translated_pr, how='previous').df,
                                snp_pr.nearest(translated_pr, how='next').df])

    adjacent_feats = adjacent_feats.groupby(['Start'])[
        ['Distance', 'Feature_name', 'Gene_name', 'Feature_description']
    ].agg(list).reset_index().T.squeeze().drop('Start')

    return adjacent_feats


def annotate_aa(snp: pd.Series, anno_pr: pr.PyRanges, anno_genome: SeqRecord) -> pd.Series:
    """For a given row in the output_df adds columns corresponding to details of the mutation"""

    aa = {}

    ref = snp['Ref']
    alt = snp['Alt']
    position = int(float(snp['Pos']))

    #     if snp['Anno_position'] == 'Undetermined':
    #         aa_columns['Mutation_type'] = 'Undetermined'
    #         aa_columns['AA_change'] = 'Undetermined'
    #         return pd.Series(aa_columns)

    if snp['Feature_type'] != 'CDS':
        aa['Mutation_type'] = np.nan
        aa['AA_change'] = np.nan

        return pd.Series(aa)

    elif len(ref) < len(alt):
        if len(ref) - len(alt) % 3 == 0:
            aa['Mutation_type'] = 'In-frame Ins'
        else:
            aa['Mutation_type'] = 'Frameshift Ins'

        aa['AA_change'] = np.nan

    elif len(ref) > len(alt):
        if len(ref) - len(alt) % 3 == 0:
            aa['Mutation_type'] = 'In-frame Del'
        else:
            aa['Mutation_type'] = 'Frameshift Del'

        aa['AA_change'] = np.nan

    else:  # len(ref) == len(alt) case
        ref_protein, mut_protein = mutate_protein(snp, anno_pr, anno_genome)

        aa['Mutation_type'], aa['AA_change'] = extract_missense(ref_protein, mut_protein)

    return pd.Series(aa)


def mutate_protein(snp: pd.Series, anno_pr: pr.PyRanges, anno_genome: SeqRecord) -> (str, str):
    """Using the annotated position, creates the WT and the mutant protein"""

    anno = anno_pr.df

    # extract snp information
    feature_name = snp['Feature_name']
    snp_coord = snp['Pos']
    ref = snp['Ref']
    alt = snp['Alt']

    # extract feature information
    feature_start = anno.loc[anno['Feature_name'] == feature_name, 'Start'].values[0]
    feature_end = anno.loc[anno['Feature_name'] == feature_name, 'End'].values[0]
    feature_strand = anno.loc[anno['Feature_name'] == feature_name, 'Strand'].values[0]

    # the genomes are zero-indexed
    gene_seq = str(anno_genome[feature_start - 1: feature_end - 1].seq)
    snp_gene_pos = snp_coord - feature_start

    ## FIX IN FUTURE:
    #     # Confirm that references bases match; with extra if statements added in for when snp occurs at end of gene
    #     if (snp_gene_pos + len(ref)) < (feature_end - feature_start):
    #         assert gene_seq[snp_gene_pos: snp_gene_pos + len(ref)] == ref
    #     elif:
    #         assert gene_seq[snp_gene_pos: snp_gene_pos]

    # generate mutant transcript
    mutated_gene_seq = gene_seq[:snp_gene_pos] + alt + gene_seq[snp_gene_pos + len(ref):]

    if feature_strand == '-':
        gene_seq = str(Seq(gene_seq).reverse_complement())
        mutated_gene_seq = str(Seq(mutated_gene_seq).reverse_complement())

    # translate transcripts
    ref_protein = SeqRecord(Seq(gene_seq)).translate().seq
    mut_protein = SeqRecord(Seq(mutated_gene_seq)).translate().seq

    return str(ref_protein), str(mut_protein)


def extract_missense(ref_protein: str, mut_protein: str) -> list:
    """For a WT and mutant protein identifies missense and silent mutations"""
    missense_aa = []

    for n in range(len(ref_protein)):
        ref_aa = ref_protein[n]
        mut_aa = mut_protein[n]
        if ref_aa != mut_aa:
            # make it 1 indexed
            missense_aa.append(f"{ref_aa}{n + 1}{mut_aa}")

    if len(missense_aa) == 0:
        return ['Silent', np.nan]

    else:
        return ['Missense', ','.join(missense_aa)]


def add_lab_ref(df, lab_snps_stringent_path, lab_snps_lenient_path):
    lab_snps_stringent = pd.read_csv(lab_snps_stringent_path)
    lab_snps_lenient = pd.read_csv(lab_snps_lenient_path)

    # make column labeling lab reference strain:
    df['Lab_reference'] = lab_snps_stringent['Strain'].values[0]

    # merge lab reference snps into snp df
    df = df.merge(lab_snps_stringent[['Ref_pos', 'Ref', 'Alt']],
                  left_on='Pos', right_on='Ref_pos', how='left', suffixes=('', '_lab_strain'))

    # generate true/false column for lab ref snps and drop all other columns
    df['Present_in_lab_reference_stringent'] = df.apply(lambda x: True if
    (x['Ref'] == x['Ref_lab_strain'] and x['Alt'] == x['Alt_lab_strain'])
    else False, axis=1)
    df = df.drop(['Ref_pos', 'Ref_lab_strain', 'Alt_lab_strain'], axis=1)

    # merge lab reference snps into snp df
    df = df.merge(lab_snps_lenient[['Ref_pos', 'Ref', 'Alt']],
                  left_on='Pos', right_on='Ref_pos', how='left', suffixes=('', '_lab_strain'))

    # generate true/false column for lab ref snps and drop all other columns
    df['Present_in_lab_reference_lenient'] = df.apply(lambda x: True if
    (x['Ref'] == x['Ref_lab_strain'] and x['Alt'] == x['Alt_lab_strain'])
    else False, axis=1)
    df = df.drop(['Ref_pos', 'Ref_lab_strain', 'Alt_lab_strain'], axis=1)

    return df


