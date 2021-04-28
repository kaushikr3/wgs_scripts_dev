import pandas as pd
import numpy as np
import os

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from Bio import GenBank

from bisect import bisect_left


def get_feature_object(coordinate, gb_features):
    """
    :param coordinate: coordinate of interest (int)
    :param gb_features: list of GenBank features
    :return: list of SeqRecord objects present in GenBank, which the passed coordinate falls inside of
    """
    feat_list = []

    for feat in gb_features:
        if coordinate in feat:  # Coordinate is simple integer
            feat_list.append(feat)

    return feat_list


def generate_annotation_series(var_coord, ref, alt, gb_feats, start_coord_dict, end_coord_dict):
    ## multiple features annotated at any given coordinate will be separated with a '/'
    feat_list = get_feature_object(var_coord, gb_feats)

    if len(feat_list) == 0:
    ## intergenic
    ## insert a check in here to manage exclusive python indexing vs. inclusive genbank notation

    else:
        # generate a row for each feature entry and then collapse if len>1
        for feat in feat_list:
            data = pd.Series(data=[feat.type, feat.qualifiers['old_locus_tag']],
                             index=['Variant_location', 'Feature_name'], name=var_coord)

            if feat.type == 'CDS':
                if 'gene' in feat.qualifiers.keys():
                    data['Variant_location'] = 'Gene'
                    data.append(pd.Series([feat.qualifiers['gene']], index=['Gene_name']))

                # evalutate insertion, deletion, frameshift insertion, frameshift deletion, missense, truncate, silent

                if len(ref) > len(alt):
                    # eval for deletion, frameshift deletion
                    aa_change = ""

                    if len(ref) - len(alt) % 3 == 0:
                        mut_type = 'In-frame Del'
                    else:
                        mut_type = 'Frameshift Del'

                elif len(ref) < len(alt):
                    aa_change = ""

                    if len(ref) - len(alt) % 3 == 0:
                        mut_type = 'In-frame Ins'
                    else:
                        mut_type = 'Frameshift Ins'

                else:
                    # eval for missense, truncate, silent
                    if 'transl_table' in feat.qualifiers.keys():
                        transl_table = int(feat.qualifiers['transl_table'])
                    else:
                        transl_table = 11

                    ref_protein = feat.qualifiers['translation']

                    ### Swap this over to generate the mutated transcript; then will compare with prot_seq
                    ## also determine aa_change here!
                    if feat.location.strand == -1:
                        ref_transcript = gb_rec.seq[
                                         int(feat.location.start):(int(feat.location.end) + 1)].reverse_complement()

                    else:
                        alt_transcript = (gb_rec.seq[int(feat.location.start): var_coord] +
                                          Seq(alt) +
                                          gb_rec.seq[var_coord + len(ref):int(feat.location.end) + 1])


        return data


def annotate_dfs(dfs, gbk_path):
    gb_rec = [rec for rec in SeqIO.parse(gbk_path, "genbank")][0]
    gb_feats = [feat for feat in gb_rec.features if ((feat.type != "gene") & (feat.type != "source"))]
    genome_seq = gb_rec.seq

    start_coord_dict = {i: int(feat.location.start) for i, feat in enumerate(gb_feats)}
    end_coord_dict = {i: int(feat.location.end) for i, feat in enumerate(gb_feats)}

    annotation_cols = ['Variant_location', 'Feature_name', 'Gene_name', 'Mutation_type', 'AA_change',
                       'Intergenic_upstream_dist', 'Intergenic_upstream_feat',
                       'Intergenic_downstream_dist', 'Intergenic_downstream_feat']
    annotated_dfs = []

    for df in dfs:
        annotate = pd.DataFrame(columns=annotation_cols)
        for i in range(len(df)):

            # function returns pd.Series named for coordinate position with index of expected annotation column names
            annotation_series = generate_annotation_series(df['POS'], df['REF'], df['ALT'],
                                                           gb_feats, genome_seq,
                                                           start_coord_dict, end_coord_dict)
            annotate.append(annotation_series[coord])

        # merge df and annotation data
        annotated_dfs.append(
            df.merge(annotate, how='left', left_on='POS', right_index=True)
        )  ## check this index/right df name situation!!!

    return annotated_dfs


### FUNCTIONS STILL TO WRITE:


def get_nearest_feature_idx(coord, start_dict, end_dict):

    return []

def format_vcf_fields(df, source):
    # returns vcf tsvs with common sense column names, and error and accuracy metric cols added in
    # any other columns worth adding??

    return df


def make_thresholded_dfs(gatk, free):
    # merges dfs, generates merged streamlined confidence metrics etc; comb through and make sure no variants
    # are duplicated due to slightly differing coordinate calls!
    df = gatk.merge(free, how="left", on="POS", suffixes=("_gatk", "_free"))

    return df


def check_lab_ref(high, low, lab_variant_csv_path):
    ### add columns for lenient and stringent lab variants

    return high, low


def write_high_conf_df(df):


## write high confidence df to .csv file


def write_low_conf_df(df):
# write variants to .xslx in series of sheets with 75% --> 99.99999% accuracy as sheet names

def main(parsed_gatk_path, parsed_freebayes_path,
         ref_gbk_path, lab_variant_csv_path):
    gatk = pd.read_csv(parsed_gatk_path, delimiter='\t')
    free = pd.read_csv(parsed_freebayes_path, delimiter='\t')

    # add accuracy and error columns, rename columns
    gatk = format_vcf_fields(gatk, source='gatk')
    free = format_vcf_fields(free, source='freebayes')

    # generate high confidence and low confidence dataframes
    high, low = make_thresholded_dfs(gatk, free)

    # adds annotation data to each df
    [high, low] = annotate_dfs([high, low], ref_gbk_path)

    # check if variant present in lab reference
    high, low = check_lab_ref(high, low, lab_variant_csv_path)

    high = high.fillna('NA')
    low = low.fillna('NA')

    # writes variant data to output files:
    write_high_conf_df(high)
    write_low_conf_df(low)