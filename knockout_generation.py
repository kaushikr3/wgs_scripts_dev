from __future__ import annotations

import argparse
import os

from Bio import SearchIO, SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter

#repl_seq_path = 'test1/KO.fa'
#repl_seq = SeqIO.read(repl_seq_path, format='fasta')
#
#
#ref_genome_path = '../../Reference/H37RvCO/H37RvCO.fasta'
#ref_genome = SeqIO.read(ref_genome_path, format='fasta')
#
def main():
    my_parser = argparse.ArgumentParser(prog='knockout_generation')
    my_parser.add_argument('-repl_seq', required=True,
                           dest='repl_seq_path', type=str)
    my_parser.add_argument('-ref_genome', required=True,
                           dest='ref_genome_path', type=str)
    my_parser.add_argument('-output_name', required=True,
                           dest='output_name', type=str)
    args = my_parser.parse_args()
    output_name = args.output_name
    repl_seq = SeqIO.read(args.repl_seq_path, format='fasta')
    create_mutant_genome(repl_seq, args.ref_genome_path, output_name)
    # Create other indexes etc for GATK
    os.system(f"bwa index {output_name}")
    os.system(f"samtools faidx {output_name}")
    os.system(
        f"java -jar ~/biotools/picard/build/libs/picard.jar CreateSequenceDictionary R={output_name} O={output_name.split('.')[0]}.dict")
    os.system(f'makeblastdb -in {output_name} -dbtype nucl')
    return None


def create_mutant_genome(repl_seq: SeqRecord, ref_genome_path: str, output_name: str):
    """With a repl_seq creates a mutant genome. Remember to remove any PmeI sites!"""
    # make zero indexed!
    up_coord = blast(repl_seq, ref_genome_path, 'upstream') - 1
    down_coord = blast(repl_seq, ref_genome_path, 'downstream') - 1
    if (up_coord == None) or (down_coord == None):
        return None
    print(up_coord, down_coord)
    ref_genome = str(SeqIO.read(ref_genome_path, 'fasta').seq)
    new_genome = ref_genome[0:up_coord] + \
        str(repl_seq.seq) + ref_genome[down_coord+1:]
    new_genome = SeqRecord(Seq(new_genome), id=output_name.split('.')[0])

    with open(output_name, "w") as o:
        SeqIO.write(new_genome, o, "fasta")
        #FastaWriter(o, wrap=0).write_file([new_genome])
    print('Sequence of new genome written to file')
    return None


def blast(repl_seq: SeqRecord, ref_genome_path: str, direction: str, query_len=200):
    if direction == 'upstream':
        query_seq = repl_seq[0:query_len]
    elif direction == 'downstream':
        query_seq = repl_seq[-query_len:]
    query_path = f"query_{direction}.fasta"
    with open(query_path, "w") as o:
        SeqIO.write(query_seq, o, "fasta")
    output_name = f"results_{direction}.xml"
    os.system(
        f"blastn -db {ref_genome_path} -query {query_path} -out {output_name} -outfmt 5")
    result_handle = open(output_name)
    blast_record = NCBIXML.read(result_handle)
    print(f"\n---{direction}---")
    print(f'#Alignments: {len(blast_record.alignments)}')
    print(f'#HSPs: {len(blast_record.alignments[0].hsps)}')
    hsp = blast_record.alignments[0].hsps[0]
    if hsp.sbjct_end < hsp.sbjct_start:
        print('The input sequence is on the complementary strand! Please provide revcomp!!')
        return None
    print(f'query_length: {query_len}')
    print(f'HSP_align_length: {hsp.align_length}')
    print(f'#Identities: {hsp.identities}')
    print(f'Query start, Query end: {hsp.query_start} , {hsp.query_end}')
    proceed = input("Proceed? Type y or n: ")
    if proceed == 'n':
        print('Aborted!')
        return None
    elif proceed == 'y':
        if direction == 'upstream':
            # 1-indexed
            return hsp.sbjct_start
        if direction == 'downstream':
            # 1-indexed
            return hsp.sbjct_end


if __name__ == "__main__":
    main()
