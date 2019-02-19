#!/usr/bin/env python3

import argparse
import pandas
import dask.dataframe as dd
from dask.multiprocessing import get


def kmerize(seq, k, stride=1, alphabet='ACGTNBDHU', complement='TGCANUHDB', delim=' '):
    """
    Convert a sequence to kmer tokens

    :param seq str: Sequence to convert
    :param k int: Size of each kmer
    :param stride int: Stride of size k window (*default: 1*)
    :param alphabet str: Alphabet of kmers (*default: 'ACGTNBDHU'*)
    :param complement str: Complement of alphabet (*default: 'TGCANUHDB'*)
    :param delim str: Delimiter of kmer token words (*default: ' '*)
    """
    dictionary = str.maketrans(alphabet, complement)
    tokens = []
    for i in range(0, len(seq) - k + 1, stride):
        kmer = seq[i:i + k].upper()
        rev_comp = kmer.translate(dictionary)[::-1]
        token = max(kmer, rev_comp)
        tokens.append(token)
    return delim.join(tokens)


def get_tokens(seqs, k, stride=1, alphabet='ACGTNBDHU', complement='TGCANUHDB', delim=' ', chunksize=100):
    """
    Convert a corpus of sequences to kmer tokens

    :param seq str: Sequence to convert
    :param k int: Size of each kmer
    :param stride int: Stride of size k window (*default: 1*)
    :param alphabet str: Alphabet of kmers (*default: 'ACGTNBDHU'*)
    :param complement str: Complement of alphabet (*default: 'TGCANUHDB'*)
    :param delim str: Delimiter of kmer token words (*default: ' '*)
    :param chunksize int: Size of chunks to partition data (*default: 100*)
    """
    dseqs = dd.from_pandas(seqs, chunksize=chunksize)
    seqs['Seq'] = dseqs.map_partitions(lambda df: df.apply(lambda row:
                                                           kmerize(row.Seq, k, stride),
                                                           axis=1)).compute(get=get)
    return seqs.rename(columns={'Seq': 'Tokens'}).drop(columns=['Strand'])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert sequence corpus to kmer token corpus.')
    parser.add_argument('input_path', type=str, help='Path to a sequence tsv file.')
    parser.add_argument('k', type=int, help='Size of each kmer.')
    parser.add_argument('--stride', dest='stride', metavar='integer', type=int, default=1,
                        help='Stride of size k window.)')
    parser.add_argument('--alphabet', dest='alphabet', metavar='string', type=str, default='ACGTNBDHU',
                        help='Alphabet of kmers.)')
    parser.add_argument('--complement', dest='complement', metavar='string', type=str, default='TGCANUHDB',
                        help='Complement of alphabet.)')
    parser.add_argument('--delim', dest='delim', metavar='string', type=str, default=' ',
                        help='Delimiter of kmer token words.)')
    parser.add_argument('--chunksize', dest='chunksize', metavar='integer', type=int, default=100,
                        help='Size of chunks to partition data.)')
    parser.add_argument('--out', dest='output_path', metavar='output_path', type=str,
                        default='corpus.tsv',
                        help='Path to corpus of kmer tokens output file (default: \'corpus.tsv\'.)')

    args = parser.parse_args()

    seqs = pandas.read_csv(args.input_path, sep='\t', dtype={'Chr': str})

    tokens = get_tokens(seqs, args.k, stride=args.stride, alphabet=args.alphabet, complement=args.complement,
                                      delim=args.delim, chunksize=args.chunksize)

    tokens.to_csv(args.output_path, sep='\t', index=False)