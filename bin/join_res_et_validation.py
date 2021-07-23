#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
import pandas as pd

# Create the parser
parser = argparse.ArgumentParser(prog='join_res_et_validation.py',
				 allow_abbrev=False,
				 usage='%(prog)s [options] <mlst_dir> <resistance_dir> <outfile>',
				 description='',
				 epilog='---')

parser.add_argument('validation_data',
		    type=str,
		    help='Place where all the results tsv:s are spread around')

parser.add_argument('nf_results',
		    type=str,
		    help='Place where all the results tsv:s are spread around')

parser.add_argument('outfile', 
		    type=str,
		    help='Place where we want our joined jsons be saved as one tsv file')


# Execute the parse_args() method
args = parser.parse_args()


val_data = Path(args.validation_data)
nf_results = Path(args.nf_results)
out = Path(args.outfile)

val_df = pd.read_csv(val_data, sep='\t')
nf_results_df = pd.read_csv(nf_results, sep='\t')

out_df = val_df.concat(right=nf_results_df, 
		       how="left", 
		       on='Prov-nummer', 
		       suffixes=["_val","_res"])

out_df.to_csv(out, index=False, sep='\t')