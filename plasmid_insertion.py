import argparse
import pandas as pd
import os

from Bio import SearchIO, SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter
from pathlib import Path


def main():
    my_parser = argparse.ArgumentParser(prog='plasmid_insertion')
    my_parser.add_argument('-plasmid_seqs', required=True,
                           dest='plasmid_seq_paths', type=str)
    my_parser.add_argument('-sites', required=True,
                           dest='sites', type=str)
    my_parser.add_argument('-strain', required=True, choices=['AL123456', 'H37RvCO', 'Erdman'],
                           dest='strain', type=str)
    my_parser.add_argument('-ref_genome', required=True,
                           dest='ref_genome_path', type=str)
    my_parser.add_argument('-output_name', required=True,
                           dest='output_name', type=str)
    my_parser.add_argument('-strain_summary_file',
                           required=True, dest='summary_file', type=str)

    args = my_parser.parse_args()
    output_name = args.output_name
    plasmid_seq_paths = args.plasmid_seq_paths.split(',')
    plasmid_seqs = [SeqIO.read(path, format='fasta')
                    for path in plasmid_seq_paths]
    sites = args.sites.split(',')

    try:
        assert len(sites) == len(plasmid_seqs)
    except AssertionError:
        print("Number of sites should  be equal to the number of plasmid sequences")
        return None

    try:
        assert set(sites).issubset(set(['tweety', 'l5', 'giles']))
    except AssertionError:
        print("Sites should be one of ['tweety', 'l5', 'giles']")
        return None

    sites_plasmids_dict = {k: v for k, v in zip(sites, plasmid_seqs)}

    # make data df
    attP_df = pd.read_csv('~/seds/wgs/wgs_scripts/metadata/attP_sites.csv')
    attB_df = pd.read_csv('~/seds/wgs/wgs_scripts/metadata/plasmid_insertion_sites.csv')
    attB_df = attB_df[(attB_df['Strain'] == args.strain)
                      & (attB_df['Site'].isin(sites))]
    data_df = pd.merge(attB_df, attP_df, how='left', on='Site')
    data_df['Plasmid_seq'] = data_df['Site'].map(sites_plasmids_dict)
    # very important sorting step, you have to insert from furthest genome position to nearest!
    data_df = data_df.sort_values(
        by="Position_in_unaltered_genome", ascending=False).reset_index()
    data_df['attP_check'] = data_df.apply(check_attP, axis=1)
    try:
        assert data_df['attP_check'].all()
    except AssertionError:
        print("Plasmid sequence(s) are incorrect. Please check if they are rotated accurately!")
        return None
    data_df['Blast_output_name'] = Path(output_name).stem + '_temp' + \
        pd.Series(data_df.index).astype(str)
    data_df['Ref_genome_path'] = args.ref_genome_path

    # make attB into seqRecord for blast
    data_df['attB_seq_start_100bp'] = data_df['attB_seq_start_100bp'].apply(
        lambda x: SeqRecord(Seq(x)))
    data_df['Blast_insertion_coord'] = data_df.apply(blast, axis=1)
    genome_to_edit = str(SeqIO.read(args.ref_genome_path, 'fasta').seq)

    # iteratively insert plasmids into genome
    for _, row in data_df.iterrows():
        new_genome = create_mutant_genome(
            row['Plasmid_seq'], genome_to_edit, row['Blast_insertion_coord'])
        genome_to_edit = new_genome

    final_genome = SeqRecord(Seq(new_genome), id=Path(output_name).stem)
    with open(output_name, "w") as o:
        SeqIO.write(final_genome, o, "fasta")
        #FastaWriter(o, wrap=0).write_file([final_genome])
    os.system(f"bwa index {output_name}")
    os.system(f"samtools faidx {output_name}")
    os.system(
        f"java -jar ~/biotools/picard/build/libs/picard.jar CreateSequenceDictionary R={output_name} O={Path(output_name).stem}.dict")
    os.system(f'makeblastdb -in {output_name} -dbtype nucl')

    # make summary file
    plasmid_names = [Path(x).stem for x in plasmid_seq_paths]
    summary_df = pd.DataFrame(
        {'Inserted_plasmid': plasmid_names, 'Plasmid_seq': plasmid_seqs, 'Site': sites})
    summary_df['Strain_name'] = Path(output_name).stem
    summary_df[['Start', 'End']] = summary_df['Plasmid_seq'].apply(
        find_plasmid_insertion_site, args=(final_genome,))
    summary_df = summary_df[['Strain_name',
                             'Inserted_plasmid', 'Site', 'Start', 'End']]
    # Since we will eventually want to check not just reads corresponding the the plasmid but also to the junction site, add 100 bp on either end
    summary_df['Start_minus100'] = summary_df['Start'] - 100
    summary_df['End_plus100'] = summary_df['End'] + 100

    summary_df = summary_df.sort_values(by='Start')
    with open(args.summary_file, 'a') as f:
        # adds header if this is the first time making the file, otherwise only appends new rows
        summary_df.to_csv(f, header=f.tell() == 0, index=False)
    return None


def check_attP(row: pd.Series):
    """Checks if the plasmid has been rotated correctly"""
    attP_seq = row['attP_sequence']
    plasmid_seq = row['Plasmid_seq']
    if str(plasmid_seq.seq)[0:len(attP_seq)] != attP_seq:
        return False
    else:
        return True


def blast(row: pd.Series):
    """Given a reference genome finds the insertion site"""
    blast_output_name = f"results_{row['Blast_output_name']}.xml"
    query_path = f"query_{row['Blast_output_name']}.fasta"
    attB_seq = row['attB_seq_start_100bp']
    with open(query_path, "w") as o:
        SeqIO.write(attB_seq, o, "fasta")
    os.system(
        f"blastn -db {row['Ref_genome_path']} -query {query_path} -out {blast_output_name} -outfmt 5")
    result_handle = open(blast_output_name)
    blast_record = NCBIXML.read(result_handle)
    try:
        assert len(blast_record.alignments) == 1
        assert len(blast_record.alignments[0].hsps) == 1
        hsp = blast_record.alignments[0].hsps[0]
        assert hsp.align_length == len(attB_seq)
        assert hsp.identities == len(attB_seq)
        # 1-indexed
        return hsp.sbjct_start
    except AssertionError:
        print("blast failed")
        return None


def create_mutant_genome(plasmid_seq: SeqRecord, ref_genome: str, blast_insertion_coord: int):
    """With a plasmid_seq creates a mutant genome"""
    # make zero indexed!
    coord = blast_insertion_coord - 1
    new_genome = ref_genome[0:coord] + \
        str(plasmid_seq.seq) + ref_genome[coord:]
    return new_genome


def find_plasmid_insertion_site(plasmid: SeqRecord, final_genome: SeqRecord):
    """After the final genome is made i.e. all plasmids are inserted, for a given plasmid finds the insertion sites"""
    final_genome = str(final_genome.seq)
    plasmid = str(plasmid.seq)
    # 0 indexed
    start = final_genome.find(plasmid)
    # make 1 indexed
    start += 1
    # end is inclusive
    end = start + len(plasmid) - 1
    return pd.Series({'Start': start, 'End': end})


if __name__ == "__main__":
    main()
