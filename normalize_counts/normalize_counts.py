#!/usr/bin/env python3

import argparse
import pandas


def normalize_length(sample_counts, gene_lengths):
	"""
	Normalize sample counts by gene length

	:param pandas.DataFrame sample_counts: Dataframe of sample feature counts
	:param pandas.Series gene_lengths: Series of gene lengths
	:return: sample feature counts normalized by length
	:rtype: pandas.DataFrame
	"""
	return sample_counts.div(gene_lengths * (10**-3), axis=0)


def normalize_reads(sample_counts):
	"""
	Normalize sample counts by reads

	:param pandas.DataFrame sample_counts: Dataframe of sample feature counts
	:return: sample feature counts normalized by number of reads
	:rtype: pandas.DataFrame
	"""
	scale_factor = sample_counts.sum(axis=0) * (10**-6)
	return sample_counts.div(scale_factor, axis=1)


def rpkm(sample_counts, gene_lengths):
	"""
	Compute RPKM (Reads Per Kilobase Million)

	:param pandas.DataFrame sample_counts: Dataframe of sample feature counts
	:param pandas.Series gene_lengths: Series of gene lengths
	:return: RPKM-normalized sample feature counts
	:rtype: pandas.DataFrame
	"""
	rpm = normalize_reads(sample_counts)
	return normalize_length(rpm, gene_lengths)


def tpm(sample_counts, gene_lengths):
	"""
	Compute TPM (Transcripts Per Million)

	:param pandas.DataFrame sample_counts: Dataframe of sample feature counts
	:param pandas.Series gene_lengths: Series of gene lengths
	:return: TPM-normalized sample feature counts
	:rtype: pandas.DataFrame
	"""
	rpk = normalize_length(sample_counts, gene_lengths)
	return normalize_reads(rpk)


function_map = {'rpkm': rpkm, 'tpm': tpm}

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Normalize feature counts.')
	parser.add_argument('counts_path', type=str, help='Path to a feature count file.')
	parser.add_argument('--out', dest='output_path', metavar='output_path', type=str,
		                default='norm_counts.tsv',
		                help='Path to normalized feature count output file (default \'norm_counts.tsv\'.)')
	parser.add_argument('--norm-type', dest='norm_type', metavar='normalization_type',
		                default=tpm, choices=function_map, type=str,
	                    help='Normalization type for feature counts (default: \'tpm\').')
	parser.add_argument('--round', dest='decimals', metavar='decimals',
		                type = int, default=3, help='Number of decimal places to round to (default: 3).')

	args = parser.parse_args()

	counts = pandas.read_csv(args.counts_path, sep='\t', header=1, index_col=0)

	length_ind = counts.columns.tolist().index('Length')
	counts.rename(columns={x: x.split('/')[-1].rstrip('.bam') for x in counts.columns[length_ind + 1:]}, inplace=True)

	norm_counts = args.norm_type(counts.iloc[:, length_ind + 1:], counts.Length)
	counts.iloc[:, length_ind + 1:] = norm_counts.round({x: args.decimals for x in norm_counts.columns})

	counts.to_csv(args.output_path, index_label='Geneid', sep='\t')