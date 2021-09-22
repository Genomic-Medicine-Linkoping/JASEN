#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Read MLST raw typing file and print out the columns of interest to STDOUT"""

import argparse
from pathlib import Path
import pandas as pd
import sys
#from io import StringIO

def print_to_stdout(*a):
	print(*a, file = sys.stdout, sep = '\t')

def print_to_stderr(*a):
	print(*a, file = sys.stderr)


def main(args):
	mlst = Path(args.mlst)
	cols = args.cols.split(',')
	(pd.read_csv(mlst, sep='\t', encoding='utf-8')[cols]
		.to_csv(sys.stdout, index=False, sep = '\t'))


if __name__ == '__main__':
	# Create the parser
	parser = argparse.ArgumentParser(prog='subset_ColsofInterest.py',
					allow_abbrev=False,
					usage='%(prog)s [options] <sample_mlst.tsv>',
					description=__doc__,
					epilog='---')

	parser.add_argument('--mlst',
				'-m',
				type=str,
				help='Path to raw mlst tsv file')

	parser.add_argument('--cols',
				'-c', 
				type=str,
				default="sample,species,sequence_type",
				help='Comma separated list of names of the columns of interest to be subseted, e.g. sample,species,sequence_type')

	# Execute the parse_args() method
	args = parser.parse_args()
	main(args)
