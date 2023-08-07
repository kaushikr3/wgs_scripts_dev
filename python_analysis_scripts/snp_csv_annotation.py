import os
import glob
import argparse
import pandas as pd
import numpy as np
import pyranges as pr

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Blast import NCBIXML

###outside of env:
# spack load /idhtq7g


# paths to references files:
if os.path.isdir('/home/nat4004'):
    genbank_dict = {
        'H37Rv': "/home/nat4004/wgs/Reference/AL123456/AL123456.gbk",
        'AL123456': "/home/nat4004/wgs/Reference/AL123456/AL123456.gbk",
    	'HN878': "/home/nat4004/wgs/Reference/NZ_CM001043/NZ_CM001043.gbk",
        'Erdman': "/home/nat4004/wgs/Reference/Erdman/Erdman.gbk",
        'Msmeg': "/home/nat4004/wgs/Reference/NC_008596/NC_008596.gb",
        'BCG': "/home/nat4004/wgs/Reference/BCG_Pasteur/BCG_Pasteur.gb"
        }
    
    
    fasta_dict = {
        'H37Rv': "/home/nat4004/wgs/Reference/AL123456/AL123456.fasta",
        'AL123456': "/home/nat4004/wgs/Reference/AL123456/AL123456.fasta",
    	'H37RvCO': "/home/nat4004/wgs/Reference/H37RvCO/H37RvCO.fasta",
    	'HN878': "/home/nat4004/wgs/Reference/NZ_CM001043/NZ_CM001043.fasta",
        'Erdman': "/home/nat4004/wgs/Reference/Erdman/Erdman.fasta",
        'Msmeg': "/home/nat4004/wgs/Reference/NC_008596/NC_008596.fasta",
        'BCG': "/home/nat4004/wgs/Reference/BCG_Pasteur/BCG_Pasteur.fa"
        }
    
    
    lab_reference_dict = {
        'North': "/home/nat4004/wgs/Reference/Lab_references/H37Rv_Jamie_North/csv/",
        'UMass': "/home/nat4004/wgs/Reference/Lab_references/H37Rv_Caro_UMass/csv/",
    	'HN878': "/home/nat4004/wgs/Reference/Lab_references/HN878_ref_Jamie/csv/",
        'Erdman': "/home/nat4004/wgs/Reference/Lab_references/Erdman_Caro/csv/",
        'BCG': "/home/nat4004/wgs/Reference/Lab_references/BCG_WT_Josh/csv/"
        }

else:
    genbank_dict = {
        'H37Rv': "~/wgs/Reference/AL123456/AL123456.gbk",
    	'HN878': "~/wgs/Reference/NZ_CM001043/NZ_CM001043.gbk",
        'Erdman': "~/wgs/Reference/Erdman/Erdman.gbk",
        'Msmeg': "~/wgs/Reference/NC_008596/NC_008596.gb",
        'BCG': "~/wgs/Reference/BCG_Pasteur/BCG_Pasteur.gb"
        }
    
    
    fasta_dict = {
        'H37Rv': "~/wgs/Reference/AL123456/AL123456.fasta",
    	'H37RvCO': "~/wgs/Reference/H37RvCO/H37RvCO.fasta",
    	'HN878': "~/wgs/Reference/NZ_CM001043/NZ_CM001043.fasta",
        'Erdman': "~/wgs/Reference/Erdman/Erdman.fasta",
        'Msmeg': "~/wgs/Reference/NC_008596/NC_008596.fasta",
        'BCG': "~/wgs/Reference/BCG_Pasteur/BCG_Pasteur.fa"
        }
    
    
    lab_reference_dict = {
        'North': "~/wgs/Reference/Lab_references/H37Rv_Jamie_North/csv/",
        'UMass': "~/wgs/Reference/Lab_references/H37Rv_Caro_UMass/csv/",
    	'HN878': "~/wgs/Reference/Lab_references/HN878_ref_Jamie/csv/",
        'Erdman': "~/wgs/Reference/Lab_references/Erdman_Caro/csv/",
        'BCG': "~/wgs/Reference/Lab_references/BCG_WT_Josh/csv/"
        }

#
# genbank_dict = {
#     'H37Rv': 'reference_files/AL123456.gbk',
#                 # 'Erdman':
#                 'Msmeg': 'reference_files/NC_008596.gb',
#                 'BCG': 'reference_files/BCG_Pasteur.gb'}
#
#
# fasta_dict = {'H37Rv': 'reference_files/H37RvCO.fasta',
#               # 'Erdman':
#               'Msmeg': 'reference_files/NC_008596.fa',
#               'BCG': 'reference_files/BCG_Pasteur.fa'}
#

# VCF Column formatting:
format_cols_dict = {
    "CHROM": "Ref_genome",
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


def main():
    """

    :return: Writes excel sheet holding annotated SNP data from parsed VCF files
    """
    parse = argparse.ArgumentParser(prog='SNP_annotation')

    parse.add_argument('-ref_strain', required=True, dest='ref_strain', type=str, 
                       choices={'H37Rv', 'Erdman', 'HN878', 'BCG', 'Msmeg', 'AL123456', 'recombinant'},
                       help='Name of reference strain used in alignemnt (H37Rv if H37RvCO used')

    parse.add_argument('-recombinant_ref', required=False, dest='recombinant_fa', type=str, nargs='?', default=False,
                       help='For recombinant alignments, pass the path to recombinant_genome.fa here')

    parse.add_argument('-recombinant_files', required=False, dest='recombinant_files', type=str, nargs='?', default=False,
                       help='For recombinant alignments, a comma seperated list of paths to any plasmids or knockout cassettes')

    parse.add_argument('--H37Rv', action='store_true', dest='h37rv_homology',
                       help='Pass to generate H37Rv homology annotation data, must also pass -blast with this arg')

    parse.add_argument('-blast', required=False, dest='blast_dir', type=str, default=None,
                       help='Path to directory to write blast files into')

    parse.add_argument('-lab_strain', required=False, dest='lab_strain', type=str, nargs='?', default=False, 
                       choices={'North', 'UMass', 'Erdman', 'HN878', 'BCG', None},
                       help='(Optional) Name of lab reference strain to add to annotation, if desired')

    parse.add_argument('-vcf', required=True, dest='vcf_dir', type=str,
                       help='Path to directory holding PARSED vcfs')

    parse.add_argument('-out', required=True, dest='out', type=str,
                       help='Path to directory to write excel files into')

    args = parse.parse_args()
    os.system("mkdir {}".format(args.out))

	# if recombinant alignment (args.refstrain = 'recombinant' and -recombinant path passed)
    if args.ref_strain == 'recombinant':
        fasta_dict.update({'recombinant': args.recombinant_fa})
        #recomb_files = args.recombinant_files.split(",")

    if args.h37rv_homology or (args.ref_strain == 'H37Rv'):
        os.system("mkdir {}".format(args.blast_dir))

    files = pd.Series([f for f in os.listdir(args.vcf_dir) if os.path.isfile(os.path.join(args.vcf_dir, f))])
    samples = files.str.split('.', expand=True)[0].str.rstrip('gatk').str.rstrip('_').unique()

    for sample in samples:
        print("Running sample: ", sample)
        vcf_dict = {}
        blast_sub_dir = None

        if args.h37rv_homology or (args.ref_strain == 'H37Rv'):
            blast_sub_dir = args.blast_dir
        #    blast_sub_dir = os.path.join(args.blast_dir, sample)
        #    os.system(f"mkdir {blast_sub_dir}")

        gatk_haploid_tsv = files[files.str.contains("{}_.+haploid".format(sample))].item()
        gatk_diploid_tsv = files[files.str.contains("{}_.+diploid".format(sample))].item()
        
        vcf_dict.update({'gatk_haploid': os.path.join(args.vcf_dir, gatk_haploid_tsv)})
        vcf_dict.update({'gatk_diploid': os.path.join(args.vcf_dir, gatk_diploid_tsv)})

        print(vcf_dict)
        df = get_full_annotated_csv(vcf_dict['gatk_haploid'], vcf_dict['gatk_diploid'],
                                    args.ref_strain, args.h37rv_homology, blast_sub_dir, args.lab_strain)

        write_excel_file(df, os.path.join(args.out, sample + "_SNP.xlsx"))


def get_full_annotated_csv(gatk_haploid_parsed_path, gatk_diploid_parsed_path,
                           reference_strain, h37rv_homology=False, blast_dir=None, lab_strain=False):
    """
    Generates fully annotated SNP df with reference annotation, H37Rv homology, and lab strain intersection

    :param gatk_haploid_parsed_path: Path to gatk haploid/stringent parsed vcf file
    :param gatk_diploid_parsed_path: Path to gatk diplid/lenient parsed vcf file
    :param reference_strain: Name of reference strain (strain aligned against -- list H37Rv NOT H37RvCO)
    :param h37rv_homology: Whether to generate and return H37Rv homology data (only used for non-H37Rv samples)
    :param blast_dir: Path to directory to store blast results in (generated for H37Rv and H37Rv homolog samples)
    :param lab_strain: Name of Lab strain to pull into annotation, if using (str)

    :return: Fully annotated SNP df, ready to be written to outfile
    """
    # open SNPs
    hgatk = pd.read_csv(gatk_haploid_parsed_path, delimiter='\t')
    dgatk = pd.read_csv(gatk_diploid_parsed_path, delimiter='\t')

    # add SNP_tool, accuracy and error columns, and rename columns
    hgatk = format_vcf_fields(hgatk, 'stringent')
    dgatk = format_vcf_fields(dgatk, 'lenient')

    # merge SNP sets with different leniency together, so subsequent operations can be performed once:
    df = pd.concat([hgatk, dgatk]).reset_index(drop=True)

    # if reference strain is recombinant, no homology info is needed, we're just looking
    # at SNPs in the cloned loci
    if reference_strain == 'recombinant':
        # later on, add in ability to limit these by coordinates of recombinant features

        df = df.fillna('NA')
        return df

    # open reference annotation data:
    ref_genbank = SeqIO.read(genbank_dict[reference_strain], "genbank")

    # if reference is H37Rv, get H37RvCO to H37Rv homology information and merge in annotation data:
    if reference_strain == 'H37Rv':
        h37rv_annotation_location = fasta_dict['H37Rv']
        h37rvCO_reference_genome = SeqIO.read(fasta_dict['H37RvCO'], format="fasta")

        df[['H37Rv_homolog_hits', 'H37Rv_homolog_%identity', 'H37Rv_homolog_position']] = df.apply(
            blast_h37rv, axis=1, args=(h37rvCO_reference_genome, blast_dir, h37rv_annotation_location))

        df = generate_annotated_df(df, "H37Rv_homolog_position", ref_genbank)


    # otherwise, annotate with given reference and generate in H37Rv homology data
    else:
        df = generate_annotated_df(df, "Pos", ref_genbank)

        if h37rv_homology:

            annotation_col_dict = {
                'Feature_type': 'H37Rv_feature_type',
                'Distance': 'H37Rv_distance_to_feature',
                'Feature_name': 'H37Rv_feature_name',
                'Gene_name': 'H37Rv_gene_name',
                'Feature_description': 'H37Rv_feature_description'
            }

            h37rv_annotation_location = fasta_dict['H37Rv']
            h37rv_gb = SeqIO.read(genbank_dict['H37Rv'], "genbank")
            h37rv_annotation = get_genbank_annotation_df(h37rv_gb)

            h37rv_feature_pr = pr.PyRanges(h37rv_annotation)
            h37rv_translated_pr = pr.PyRanges(h37rv_annotation[h37rv_annotation['Feature_type'] != 'regulatory'][
                                                  h37rv_annotation['Feature_type'] != 'misc_feature'])

            df[['H37Rv_homolog_hits', 'H37Rv_homolog_%identity', 'H37Rv_homolog_position']] = df.apply(
                blast_h37rv, axis=1, args=(ref_genbank, blast_dir, h37rv_annotation_location))

            position_col_name = "H37Rv_homolog_position"
            h37rv_homolog_columns = list(annotation_col_dict.values())

            # pull in snp annotation data:
            df[h37rv_homolog_columns] = df.apply(
                lambda x: annotate_snp(x, position_col_name, h37rv_feature_pr, h37rv_translated_pr,
                                       annotation_col_dict), axis=1)

            # convert list entries to strings
            df[h37rv_homolog_columns] = df[h37rv_homolog_columns].applymap(
                lambda x: ";".join(map(str, x)) if isinstance(x, list) else x)

            # simplify feature_type column; drop {};{} field unless there are actually two different feature types
            df.loc[df['H37Rv_feature_type'].str.contains(';'), 'H37Rv_feature_type'] = df.loc[
                df['H37Rv_feature_type'].str.contains(';'), 'H37Rv_feature_type'].str.split(';', expand=True).apply(
                lambda x: x[0] if (x[0] == x[1] or pd.isna(x[1]))
                else "{};{}".format(x[0], x[1]), axis=1)

    # # add lab_reference_info
    if lab_strain:
        df = add_lab_ref(df, lab_reference_dict[lab_strain])

    df = df.fillna('NA')

    return df


def generate_annotated_df(df, position_col_name, ref_genbank):
    """

    :param df: df with all SNPs
    :param position_col_name: name of Position column to use when looking up annotation data
    :param ref_genbank: reference genbank file
    :return: DF with in depth annotation info for each SNP
    """

    # generate annotation df and pr objects:
    annotation = get_genbank_annotation_df(ref_genbank)
    anno_pr = pr.PyRanges(annotation)
    translated_pr = pr.PyRanges(annotation[annotation['Feature_type'] != 'regulatory'][
                                    annotation['Feature_type'] != 'misc_feature'])

    # pull in snp annotation data:
    df[['Feature_type', 'Distance_to_feature', 'Feature_name', 'Gene_name', 'Feature_description']] = df.apply(
        lambda x: annotate_snp(x, position_col_name, anno_pr, translated_pr), axis=1)

    # separate out data so that any snp with multiple features can be properly annotated
    cols_to_explode = ['Feature_type', 'Distance_to_feature', 'Feature_name', 'Gene_name', 'Feature_description']
    df = pd.DataFrame(df.set_index([c for c in df.columns if c not in cols_to_explode])
                      ).explode(cols_to_explode).reset_index()

    # annotate mutation types in coding sequences
    df[['Mutation_type', 'AA_change']] = df.apply(
        lambda x: annotate_aa(x, position_col_name, anno_pr, ref_genbank), axis=1)

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


def get_error_probability(quality):
    """
    :param quality: phred-scaled estimate of confidence that a polymorphism exists at a given genome position
    :return: probability that the (variant) base has been called incorrectly by the sequencer
    """
    error = 10 ** (-quality / 10)

    return error


def format_vcf_fields(df, stringency):

    column_check = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'sample.GT', 'sample.GQ', 'sample.PL', 'sample.AD',
                    'sample.DP']
    #assert df.columns.tolist() == column_check  # vcf df in unexpected format

    # split allele depth column into reference allele count and alternate allele count columns
    df[["Ref_reads", "Alt_reads"]] = df['sample.AD'].str.split(',', expand=True)[[0, 1]]
    df = df.drop(columns=['sample.PL', 'sample.AD'])

    # dropping this in here because excel files are being written with these cols as strings
    df['Alt_reads'] = df['Alt_reads'].astype(int)
    df['Ref_reads'] = df['Ref_reads'].astype(int)

    df['%Reads_alt_allele'] = df['Alt_reads'] / df['sample.DP']
    df['%Reads_alt_allele'] = df['%Reads_alt_allele'].round(4)

    df = df.rename(columns=format_cols_dict)

    df['Error_prob'] = df['Qual'].apply(lambda x: get_error_probability(x))
    df['Accuracy_prob'] = 1 - df['Error_prob']

    # add in column labeling if SNPs are stringent or lenient:
    df['Stringency'] = stringency

    # reorder columns and drop composite_hits:
    df = df[['Ref_genome', 'Pos', 'Ref', 'Alt', 'Qual', 'Error_prob', 'Accuracy_prob',
             'GT', 'GQ', 'DP', 'Ref_reads', 'Alt_reads', '%Reads_alt_allele', 'Stringency']]

    # returns vcf tsvs with common sense column names, and error and accuracy metric cols added in
    # any other columns worth adding??

    return df


def get_genbank_annotation_df(genome: SeqRecord) -> pd.DataFrame:
    """
    :param genome: SeqRecord object from reference genbank file
    :return: df with features, feature descriptors and coordinates from genbank file
    """
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
    else:
        feats['Feature_name'] = feats['Feature_tag']

    # make this 1-indexed; biopython reads in genbanks 0-based w/ non-inclusvie end coord
    feats['Start'] = feats['Start'] + 1
    feats['End'] = feats['End']

    return feats


def annotate_snp(snp, position_col_name, anno_pr, translated_pr, annotation_col_dict=None):
    """

    :param snp: SNP row being annotated (pd.Series)
    :param position_col_name: Name of column to use as reference coordinate (either calls VCF Ref value,
            or H37Rv homolog coordinated) (str)
    :param anno_pr: pr.PyRanges object holding all annotation data in annotation (pr.PyRanges)
    :param translated_pr: pr.PyRanges object holding all annotation data for translated loci in annotation (pr.PyRanges)
    :param annotation_col_dict: (optional) Dictionary mappinng default column names to alternate
            (H37Rv homolog) column names (dict)
    :return: pd.Series with full annotation information for a SNP (holds Feature_type, Distance, Feature_name,
            Gene_name, and Feature Description) (pd.Series)
    """
    snp_coord = snp[position_col_name]
    if snp_coord == 'Undetermined':
        snp_annotation = pd.Series({k: v for k, v in
                                    zip(['Feature_type', 'Distance', 'Feature_name', 'Gene_name',
                                         'Feature_description'],
                                        ['Undetermined'] * 5)})

        if annotation_col_dict:
            return snp_annotation.rename(index=annotation_col_dict)
        else:
            return snp_annotation

    snp_pr = pr.PyRanges(chromosomes=anno_pr.chromosomes[0], starts=[snp_coord], ends=[snp_coord])
    intersect = anno_pr.intersect(snp_pr).df

    # return flanking features if SNP is intergenic or occurs in a regulatory feature
    if len(intersect) == 0:
        snp_annotation = get_adjacent_features(anno_pr, snp_pr)
        snp_annotation = pd.concat([pd.Series({'Feature_type': len(snp_annotation['Distance'])*['Intergenic']}),
                                    snp_annotation])

    # if only intersecting features are regulatory or misc_features, return nearest features instead
    elif set(intersect['Feature_type']).union({'regulatory', 'misc_feature'}) == {'misc_feature', 'regulatory'}:
        snp_annotation = get_adjacent_features(translated_pr, snp_pr)
        snp_annotation = pd.concat(
            [pd.Series({'Feature_type': 
                       ["Intergenic ({})".format(",".join(intersect['Feature_description'])) for i in range(len(
                                 snp_annotation['Distance']))]}),
             snp_annotation])

    # if one overlapping feature return that feature's details
    elif len(intersect) == 1:
        intersect['Distance'] = np.nan
        snp_annotation = intersect[
            ['Feature_type', 'Distance', 'Feature_name', 'Gene_name', 'Feature_description']].T.squeeze()

    # else more than one feature; squeeze into one row
    else:
        intersect['Distance'] = np.nan
        snp_annotation = intersect.groupby(['Start'])[
            ['Feature_type', 'Distance', 'Feature_name', 'Gene_name', 'Feature_description']
        ].agg(list).reset_index().T.squeeze().drop('Start')

    if annotation_col_dict:
        return snp_annotation.rename(index=annotation_col_dict)
    else:
        return snp_annotation


def get_adjacent_features(anno_pr: pr.PyRanges, snp_pr: pr.PyRanges) -> pd.Series:
    """
    Extracts information of neighbouring features given an intergenic genomic position

    :param anno_pr: pr.PyRanges object holding annotation data from annotation genbank file; expected to be all
                annotation data, or annotations for all translated loci in reference
    :param snp_pr: pr.PyRanges object holding location of SNP
    :return: pd.Series holding Distance from SNP to annotated features, as well as those features' names, gene_names,
                and descriptions (pd.Series)
    """

    # pull annotation data on previous and next features:
    adjacent_feats = pd.concat([snp_pr.nearest(anno_pr, how='previous').df,
                                snp_pr.nearest(anno_pr, how='next').df])

    adjacent_feats = adjacent_feats.groupby(['Start'])[
        ['Distance', 'Feature_name', 'Gene_name', 'Feature_description']
    ].agg(list).reset_index().T.squeeze().drop('Start')

    return adjacent_feats


def blast_h37rv(row: pd.Series, ref_genome: SeqRecord, blast_dir_name: str,
                annotation_ref_location: str, query_len=400) -> pd.Series:
    """
    For each row, uses the reference position to create a query sequence of length query length.
    This is then blast-ed against the annotated genome. The alignment is verified and details are returned as a Series
    :param row: SNP row entry
    :param ref_genome: SeqRecord object of reference genome (this is the sequence reads were aligned against)
    :param blast_dir_name: Name of subdirectory to store blast files in
    :param annotation_ref_location: path to Blast DB for annotation strain
    :param query_len: length of query to blast
    :return: pd.Series - (#hits, %identity of best hit, and annotation position homologous to SNP position)
    """

    # Note that reference genome is zero-indexed and ref position is 1-indexed
    ref_position = row['Pos']
    # ref = row['Ref']

    file_name = str(ref_position)
    query_len = int(query_len)
    half_query_len = query_len // 2

    query_seq = ref_genome[ref_position - 1 - half_query_len: ref_position + half_query_len]

    query_path = os.path.join(blast_dir_name, "temp_query.fasta")
 
    with open(query_path, "w") as o:
        SeqIO.write(query_seq, o, "fasta")

    output_path = os.path.join(blast_dir_name, "temp_results.xml")
    #output_path = os.path.join(blast_dir_name, "results_{}.xml".format(file_name))

    os.system("blastn -db {ref} -query {query} -out {out} -outfmt 5".format(
            ref=annotation_ref_location, query=query_path, out=output_path))

    result_handle = open(output_path)
    blast_record = NCBIXML.read(result_handle)
    col_names = ['#HSPs', 'Percent_identity', 'Anno_position']
    n_alignments, identities, anno_position = verify_homolog_alignment(
        blast_record, query_len)

    return pd.Series({k: v for k, v in zip(col_names, [n_alignments, identities, anno_position])})


def verify_homolog_alignment(blast_record, query_len: int) -> (int, int, int):
    """

    :param blast_record: XML file object that's been read in, holding information about results of blasting some query
                        flanking a SNP against the annotation genome of interest
    :param query_len: Length of query blasted against ref (int)
    :return: Tuple holding number of hits in blast record, %identity of best hit, and coordinate in hit corresponding to
            SNP coordinate (int, int, int)
    """
    n_hsps = 0
    try:
        assert len(blast_record.alignments) > 0

        # Holdover from Anisha's code -- not sure if needed
        if len(blast_record.alignments) > 1:
            print('MORE THAN ONE ALIGNMENT!; BLAST ARG ERROR')
            raise ValueError

        n_hsps = len(blast_record.alignments[0].hsps)
        hsp = blast_record.alignments[0].hsps[0]

        assert hsp.align_length > 0.9 * query_len
        identities = hsp.identities * 100 / query_len
        # I think remove the +1 here because it accounts for exlusive outer indexing, but we're not indexing a range
        #return n_hsps, identities, hsp.sbjct_end + 1 - (query_len // 2)
        return n_hsps, identities, hsp.sbjct_end - (query_len // 2)
    except AssertionError:
        return n_hsps, 'Undetermined', 'Undetermined'


def annotate_aa(snp: pd.Series, position_col_name: str, anno_pr: pr.PyRanges, anno_genome: SeqRecord) -> pd.Series:
    """

    :param snp: Row with SNP info from SNP df (pd.Series)
    :param position_col_name: Name of column to use for pulling reference information
                                (Either ref genome or H37Rv genome coordinate) (str)
    :param anno_pr: pr.PyRanges object holding all annotation information
    :param anno_genome: annotation genome genbank file (SeqRecord)
    :return: Series holding information on Mutation type and specific AA change
    """
    """For a given row in the output_df adds columns corresponding to details of the mutation, 
    this fuction assumes position is passed 1-indexed, and switches to to 0-index for calculations 
    """

    aa = {}
    mut_type = []
    aa_change = []

    if snp[position_col_name] == 'Undetermined':
        aa['Mutation_type'] = 'Undetermined'
        aa['AA_change'] = 'Undetermined'
        return pd.Series(aa)

    ref = snp['Ref']
    alternates = snp['Alt'].split(",")
    snp_coord = int(float(snp[position_col_name])) - 1

    if snp['Feature_type'] != 'CDS':
        aa['Mutation_type'] = np.nan
        aa['AA_change'] = np.nan

        return pd.Series(aa)

    for alt in alternates:
        alt = alt.replace("*", "")

        if len(ref) < len(alt):
            if len(ref) - len(alt) % 3 == 0:
                mut_type.append('In-frame Ins')
            else:
                mut_type.append('Frameshift Ins')

            aa_change.append("NA")

        elif len(ref) > len(alt):
            if len(ref) - len(alt) % 3 == 0:
                mut_type.append('In-frame Del')
            else:
                mut_type.append('Frameshift Del')

            aa_change.append("NA")

        else:  # len(ref) == len(alt) case
            ref_protein, mut_protein = mutate_protein(snp_coord, ref, alt, snp['Feature_name'], anno_pr, anno_genome)
            alt_mutation_type, alt_aa_change = extract_missense(ref_protein, mut_protein)

            mut_type.append(alt_mutation_type)
            aa_change.append(alt_aa_change)

    aa.update({'Mutation_type': ",".join(mut_type), 'AA_change': ",".join(aa_change)})
    return pd.Series(aa)


def mutate_protein(snp_coord: int, ref: str, alt: str, feature_name: str, 
                   anno_pr: pr.PyRanges, anno_genome: SeqRecord) -> (str, str):
    """

    :param snp: Row from SNP df (pd.Series)
    :param snp_coord: SNP coordinate -- may be SNP position, or H37Rv homolog position (int)
    :param alt: alternate bases (str)
    :param anno_pr: annotation information (pr.PyRanges object)
    :param anno_genome: reference genome SeqRecord from Genbank file
    :return: (reference protein sequence, mutant protein sequence) (str, str)
    """
    """Using the annotated position, creates the WT and the mutant protein
    assumes position value passed is 0-inded 
    """

    anno = anno_pr.df

    # extract snp information
    #feature_name = snp['Feature_name']
    #ref = snp['Ref']
    # snp_coord = snp['Pos']
    # alt = snp['Alt']

    # extract feature information
	# swapping for 0-indexed non-inclusive indexing:
    feature_start = anno.loc[anno['Feature_name'] == feature_name, 'Start'].values[0] - 1
    feature_end = anno.loc[anno['Feature_name'] == feature_name, 'End'].values[0]
    feature_strand = anno.loc[anno['Feature_name'] == feature_name, 'Strand'].values[0]

    # the genomes are zero-indexed
    gene_seq = str(anno_genome[feature_start: feature_end].seq)
    snp_gene_pos = snp_coord - feature_start

    # FIX IN FUTURE:
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
    ref_protein = SeqRecord(Seq(gene_seq)).translate(table=11).seq
    mut_protein = SeqRecord(Seq(mutated_gene_seq)).translate(table=11).seq

    return str(ref_protein), str(mut_protein)


def extract_missense(ref_protein: str, mut_protein: str) -> list:
    """

    :param ref_protein: translated reference (WT) protein sequence (str)
    :param mut_protein: translated mutant protein sequence, accounting for SNP (str)
    :return: Type of mutation (str)
    """
    missense_aa = []

    for n in range(len(ref_protein)):
        ref_aa = ref_protein[n]
        mut_aa = mut_protein[n]
        if ref_aa != mut_aa:
            # make it 1 indexed
            missense_aa.append("{ref_aa}{coord}{mut_aa}".format(ref_aa=ref_aa, coord=(n+1), mut_aa=mut_aa))

    if len(missense_aa) == 0:
        return ['Silent', 'NA']

    else:
        return ['Missense', ','.join(missense_aa)]


def add_lab_ref(df, lab_variant_csv_path):
    """

    :param df: SNP dataframe (pd.DataFrame)
    :param lab_variant_csv_path: path to directory hold lab reference SNP csvs (str)
    :return: df with new columns identifying SNPs presence in Lab reference strains merged in (pd.DataFrame)
    """
    lab_snps_stringent = pd.read_csv(
        glob.glob(os.path.join(lab_variant_csv_path, '*haploid.stringent.annotated.csv'))[0])
    lab_snps_lenient = pd.read_csv(
        glob.glob(os.path.join(lab_variant_csv_path, '*diploid.lenient.annotated.csv'))[0])

    # make column labeling lab reference strain:
    df['Lab_reference'] = lab_snps_stringent['Strain'].values[0]

    # merge stringent lab reference snps into snp df
    df = df.merge(lab_snps_stringent[['Ref_pos', 'Ref', 'Alt']],
                  left_on='Pos', right_on='Ref_pos', how='left', suffixes=('', '_lab_strain'))

    # generate true/false column for lab ref snps and drop all other columns
    df['Present_in_lab_reference_stringent'] = df.apply(
        lambda x: True if (x['Ref'] == x['Ref_lab_strain'] and x['Alt'] == x['Alt_lab_strain']) else False, axis=1)
    df = df.drop(['Ref_pos', 'Ref_lab_strain', 'Alt_lab_strain'], axis=1)

    # merge lenient lab reference snps into snp df
    df = df.merge(lab_snps_lenient[['Ref_pos', 'Ref', 'Alt']],
                  left_on='Pos', right_on='Ref_pos', how='left', suffixes=('', '_lab_strain'))

    # generate true/false column for lab ref snps and drop all other columns
    df['Present_in_lab_reference_lenient'] = df.apply(
        lambda x: True if (x['Ref'] == x['Ref_lab_strain'] and x['Alt'] == x['Alt_lab_strain']) else False, axis=1)
    df = df.drop(['Ref_pos', 'Ref_lab_strain', 'Alt_lab_strain'], axis=1)

    return df


def write_excel_file(df, filename):
    """
    Writes SNPs to Excel file with 5 sheets holding 5 levels of confidence that SNP was called correctly
    :param df: fully annotated SNP df ready to be written to Excel file (pd.DataFrame)
    :param filename: filename dest (str)
    :return: None
    """

    dfs_to_write = {
        'Stringent': df[df['Stringency'] == 'stringent'].drop(columns='Stringency').sort_values(by='Pos'),
        'Lenient': df[df['Stringency'] == 'lenient'].drop(columns='Stringency').sort_values(by='Pos')
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

