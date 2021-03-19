import pandas as pd
import pyranges as pr
from Bio import SearchIO, SeqIO
from Bio.SeqRecord import SeqRecord


def build_features(genome: SeqRecord) -> pd.DataFrame:
    """Extracts relevant features from the genome with location"""
    selected_features = []
    for feature in genome.features:
        feature_dict = {'Chromosome': genome.name}
        if feature.type == 'CDS':
            feature_dict['Feature_name'] = feature.qualifiers['locus_tag'][0]
        # for inserted plasmids etc
        if feature.type == "misc_feature":
            feature_dict['Feature_name'] = feature.qualifiers['locus_tag'][0]
        if feature.type == 'tRNA':
            feature_dict['Feature_name'] = feature.qualifiers['product'][0]
        if (feature.type == 'ncRNA') or (feature.type == 'misc_RNA'):
            feature_dict['Feature_name'] = feature.qualifiers['gene'][0]
        feature_dict['Feature_type'] = feature.type
        feature_dict['Start'] = int(feature.location.start)
        feature_dict['End'] = int(feature.location.end)
        if feature.strand == 1:
            feature_dict['Strand'] = '+'
        elif feature.strand == -1:
            feature_dict['Strand'] = '-'
        try:
            feature_dict['Feature_name']
            selected_features.append(feature_dict)
        except KeyError:
            pass
    # make this 1-indexed
    out = pd.DataFrame(selected_features)
    out['Start'] = out['Start']+1
    out['End'] = out['End']+1
    return out


anno_genome = SeqIO.read(
    './metadata/genbank_files/NZCM001043.gbk', "genbank")

selected_features = build_features(anno_genome)

print(selected_features['Feature_type'].value_counts())

selected_features.to_csv(
    './metadata/' + anno_genome.name + '_selected_features.csv', index=False)
