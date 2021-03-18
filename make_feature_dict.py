import csv
import pandas as pd
import pyranges as pr
from Bio import SearchIO, SeqIO
from Bio.SeqRecord import SeqRecord

# anno=SeqIO.read('../ref_genomes/BCG-Pasteur/BCG_Pasteur.gb',"genbank")
#anno=SeqIO.read("metadata/genbank_files/NC_008769.1.gb", "genbank")
anno=SeqIO.read("metadata/genbank_files/NZCM001043.gbk", "genbank")

"""
WANT TO GENERATE A CSV CONTANING:
    GENE ID, GENE NAME, SEQ TYPE
`
"""

with open('NZ_CM001043_gene_details.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["ORF,Feature_gene_name,Feature_description"])

    for feat in anno.features:

        if feat.type == "gene":
            gene_id=feat.qualifiers['locus_tag'][0]
            gene_desc="gene"

            if 'gene' in feat.qualifiers.keys():
                gene_name=feat.qualifiers['gene'][0]
            else:
                gene_name=feat.qualifiers["locus_tag"][0]

            writer.writerow([gene_id, gene_name, gene_desc])

        if feat.type == "CDS":
            gene_id=feat.qualifiers['locus_tag'][0]
            gene_name=feat.qualifiers['locus_tag'][0]
            gene_desc="CDS"

            writer.writerow([gene_id, gene_name, gene_desc])
