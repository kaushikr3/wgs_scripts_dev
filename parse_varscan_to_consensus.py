import pandas as pd
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse


def parse_varscan_file(varscan_filepath):
    df = pd.read_csv(varscan_filepath, delimiter='\t')
    df = df[['Ref', 'Var', 'Freq']]
    df['Freq'] = df['Freq'].str.rstrip('%')
    df['Freq'] = df['Freq'].replace('-', 0).astype(float)

    for i in df.index:
        if df.loc[i, 'Freq'] > 50:
            var_base = df.loc[i, 'Var']
            if len(var_base) == 1:
                df.loc[i,'Ref'] = df.loc[i, 'Var']

    seq = ''.join(df['Ref'].tolist())
    print("Seq length: ", len(seq))
    return seq


def write_fasta(seq, output):
    fasta_name = os.path.split(output)[-1].rstrip('fasta').rstrip('.')
    out_seq = SeqRecord(Seq(seq), id=fasta_name)
    with open(output, "w") as o:
        SeqIO.write(out_seq, o, "fasta")
    print("{} written".format(output))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', dest='varscan_file', type=str, required=True,
                        help='pre-edited (:/\t) varscan file')
    parser.add_argument('--out', dest='fasta_name', type=str, required=True,
                        help='fasta output name stem')

    args = parser.parse_args()
    seq = parse_varscan_file(args.varscan_file)
    write_fasta(seq, args.fasta_name)


if __name__ == '__main__':
    main()
