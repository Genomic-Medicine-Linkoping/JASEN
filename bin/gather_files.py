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
parser.add_argument('--pattern',
		    '-p',
		    type=str,
		    default='**/motif_report*tsv',
		    help='Recursive globbing pattern to search files with')

parser.add_argument('--suffix',
		    '-s',
		    type=str,
		    default='.tsv',
		    help='File name suffix to use for when saving the copied files')

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

# If output dir is missing create it
if not gathered_tsvs.is_dir():
	gathered_tsvs.mkdir(parents=True, exist_ok=True)

for path in sorted(all_results.rglob(args.pattern)):
	id = path.resolve().parent.parent.name
	motif_report = path.resolve()
	stem = path.stem
	new_filename = id + "_" + stem + args.suffix
	out_file_path = gathered_tsvs.resolve()/new_filename
	copy(motif_report, out_file_path)
