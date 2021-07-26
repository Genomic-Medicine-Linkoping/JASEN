#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
import pandas as pd

# Create the parser
parser = argparse.ArgumentParser(prog='compare_resistance.py',
				    allow_abbrev=False,
				    usage='%(prog)s [options] <validataion.tsv> <JASEN_results.tsv> <outfile.tsv>',
				    description='Join tables and compare resistance genes',
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


# Path to where the results are
validation = Path(args.validation)
# Path where we want our results to end up
results = Path(args.results)
# Path to the joined tsv output file 
outfile = Path(args.outfile)

val = pd.read_csv(validation, sep='\t', encoding='utf-16', encoding_errors='ignore')
res = pd.read_csv(results, sep='\t', encoding='utf-8', encoding_errors='ignore')

# print(val.columns)

merged = val[['Prov-nummer', 'Art', 'Kända resistensgener samt PVL']].merge(right=res[["Prov-nummer", "Kända resistensgener samt PVL"]], 
									      how='left', 
									      on='Prov-nummer', 
									      suffixes=[" (validering)", " (resultat)"])

merged.to_csv(outfile, sep='\t', index=False)
