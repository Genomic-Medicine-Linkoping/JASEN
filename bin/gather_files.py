#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""This script gathers files matched with <pattern> recursively from indir and copies them with <suffix> to outdir."""

import argparse
from pathlib import Path
from shutil import copy

def filenamer(id, stem, suffix, includeID):
	"""Form a filename string from id (if includeID == yes) stem and suffix"""
	if(includeID == 'y'):
		return id + "_" + stem + suffix
	elif(includeID == 'n'):
		return stem + suffix
	else:
		assert False, "Alternatives are only 'y' or 'n'."

def make_dir(dir):
	"""Create a directory if it doesn't yet exist"""
	if not dir.is_dir():
		dir.mkdir(parents=True, exist_ok=True)



def main(args):
	# Path to where the results are
	all_results = Path(args.indir)
	# Path where we want our results to end up
	gathered_tsvs = Path(args.outdir)
	make_dir(gathered_tsvs)
	# Find spread out results files and gather them to one output dir
	for path in sorted(all_results.rglob(args.pattern)):
		id = Path()
		if (args.ariba_mlst == "y"):
			id = path.resolve().parent.parent.parent.name
		else:
			id = path.resolve().parent.parent.name
		current_file_path = path.resolve()
		stem = path.stem
		new_filename = filenamer(id, stem, args.suffix, args.includeID)
		out_file_path = gathered_tsvs.resolve()/new_filename
		copy(current_file_path, out_file_path)

if __name__ == '__main__':
	# Create the parser
	parser = argparse.ArgumentParser(prog='gather_files.py',
					allow_abbrev=False,
					usage='%(prog)s [options] <indir> <outdir>',
					description='Gather files of interest from <indir> to <outdir> using "pattern" and saving with "suffix"',
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

	parser.add_argument('outdir', 	# Path to where the results are
			type=str,
			help='Where the gathered files should be copied.')

	parser.add_argument('--includeID',
			    '-i',
			    type=str,
				choices=['y','n'],
			    default='n', # When searching ariba mlst results, use yes
			    help='Whether to include ID part in the file name or not')

	parser.add_argument('--ariba-mlst',
			    '-a',
			    type=str,
				choices=['y','n'],
			    default="n", # Use 'y' with ariba mlst tsv files
			    help='Are you searching for ariba mlst files from JASEN pipeline results output dir?')


	# Execute the parse_args() method
	args = parser.parse_args()
	main(args)