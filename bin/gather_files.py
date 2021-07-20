#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
from shutil import copy

# Create the parser
parser = argparse.ArgumentParser(prog='gather_files.py',
				    allow_abbrev=False,
				    usage='%(prog)s [options] <indir> <outdir>',
				    description='List the content of a folder',
				    epilog='---')

parser.add_argument('indir',
		    type=str,
		    help='Place where all the results tsv:s are spread around')

parser.add_argument('outdir', 
		    type=str,
		    help='Place where we want our renamed tsv:s')

# Execute the parse_args() method
args = parser.parse_args()

# Path to where the results are
all_results = Path(args.indir)
# Path where we want our results to end up
gathered_tsvs = Path(args.outdir)


for path in sorted(all_results.rglob('**/motif_report*tsv')):
	id = path.resolve().parent.parent.name
	motif_report = path.resolve()
	stem = path.stem
	new_filename = id + "_" + stem + ".tsv"
	copy(motif_report, gathered_tsvs/new_filename)
