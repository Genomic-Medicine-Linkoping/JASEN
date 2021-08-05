#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Join tables and compare MLST typing"""

import argparse
from pathlib import Path
import pandas as pd

def main(args):
	# Path to where the results are
	validation = Path(args.validation)
	# Path where we want our results to end up
	results = Path(args.results)
	# Path to the joined tsv output file 
	outfile = Path(args.outfile)

	val = pd.read_csv(validation, sep='\t', encoding='utf-16', encoding_errors='ignore')
	res = pd.read_csv(results, sep='\t', encoding='utf-8', encoding_errors='ignore')

	merged = val[['Prov-nummer', 'Art', 'Känd MLST eller Spa-typ']].merge(right=res[["Prov-nummer", "Känd MLST eller Spa-typ"]], 
											how='left', 
											on='Prov-nummer', 
											suffixes=[" (validering)", " (resultat)"])

	merged.to_csv(outfile, sep='\t', index=False)


if __name__ == '__main__':
	# Create the parser
	parser = argparse.ArgumentParser(prog='compare_mlst.py',
					allow_abbrev=False,
					usage='%(prog)s [options] <validation.tsv> <JASEN_results.tsv>',
					description=__doc__,
					epilog='---')

	parser.add_argument('validation',
				type=str,
				help='Path to validation tsv file')

	parser.add_argument('results', 
				type=str,
				help='Path to JASEN results tsv file')

	parser.add_argument('outfile', 
				type=str,
				help='Path to joined output tsv file')

	# Execute the parse_args() method
	args = parser.parse_args()
	main(args)