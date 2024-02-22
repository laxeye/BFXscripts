#!/usr/bin/env python3
import argparse
import glob
import os.path
import re
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord


def parse_args():
	parser = argparse.ArgumentParser(description='Combine alignments.')

	parser.add_argument('-i', '--input', help='Input folder')
	parser.add_argument('-o', '--output', help='Output FASTA')
	'''parser.add_argument('-', '--', help='')'''

	args = parser.parse_args()

	assert os.path.isdir(args.input), 'Not a valid directory'
	args.input = os.path.abspath(args.input)

	return args


def combine_alignments(input_dir, out_fasta):
	aligns_d = dict()
	ncbi_patt = re.compile(r'GC[AF]_\d+\.\d')
	for file in glob.glob(os.path.join(input_dir,'*.f*a')):
		try:
			align = AlignIO.read(file, "fasta")
			for seq in align:
				#marker, genome = seq.id.split('.', 1)
				m = ncbi_patt.search(seq.id)
				if m:
					genome = m.group()
				else:
					genome = seq.id
				aligns_d[genome] = aligns_d.get(genome, '') + seq.seq
		except:
			print(f'Empty or invalid file: {file}')
	with open(out_fasta, 'w') as out_h:
		for k,v in aligns_d.items():
			print(f'>{k}\n{v}', file=out_h)


def main():
	args = parse_args()
	combine_alignments(args.input, args.output)


if __name__ == '__main__':
	main()
