#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""This script gathers tsv files from indir and joins them and outputs them to outdir."""

import argparse
from pathlib import Path
from shutil import copy
import pandas as pd

def main(args):
	# Path to where the tables are
	all_tables = Path(args.indir)
	# Path where we want our joined table to end up
	outfile = Path(args.outfile)
	tables_list = []
	# Loop through the table files
	for path in sorted(all_tables.glob(args.pattern)):
		tbl = pd.read_csv(path, sep=args.separator)
		
		# Parse out species and sample names
		fn_parts = str(path.stem).split('_')
		tbl['Sample'] = fn_parts[3]
		tbl['Species'] = fn_parts[1] + ' ' + fn_parts[2]
		cols = tbl.columns.tolist()
		# Move last two elements to first
		cols = cols[-2:] + cols[:-2]
		# Move last element to third place
		last = cols[-1]
		cols.insert(2, last)
		cols = cols[:-1]
		tbl = tbl[cols]
		# Add table to the list
		tables_list.append(tbl)
	pd.concat(tables_list, ignore_index=True).to_csv(outfile, sep=args.separator, index=False)

if __name__ == '__main__':
	# Create the parser
	parser = argparse.ArgumentParser(prog='join_tables.py',
					 allow_abbrev=False,
					 usage='%(prog)s [options] <indir> <outfile>',
					 description=__doc__,
					 epilog='---')

	parser.add_argument('--pattern',
			    '-p',
			    type=str,
			    default='*.tsv',
			    help='Globbing pattern to search files with')

	parser.add_argument('--separator',
			    '-s',
			    type=str,
			    default='\t',
			    help='Delimiter in the tables')	

	parser.add_argument('indir',
			    type=str,
			    help='Place where all the results tsv:s are')

	parser.add_argument('outfile',
			    type=str,
			    help='Where the joined table should be copied.')


	# Execute the parse_args() method
	args = parser.parse_args()
	main(args)