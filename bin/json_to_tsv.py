#!/usr/bin/env python

'''
Copyright 2021 Isak Sylvin
This source code file is part of  genomic-medicine-sweden/JASEN github repository 
and is licensed under the GPL-3.0. 
For the full license text, please see the file LICENSE 
or visit https://www.gnu.org/licenses/gpl-3.0.en.html
16 Sep 2021 - Modified by Lauri Mesilaakso.
'''

'''
    Usage:   python json_to_tsv.py <infile> <outfile>
    Example: python json_to_tsv.py transposed_report.json report.tsv
'''

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='This script converts a \
    json file (.json) into tab-separated report format (.tsv).')

parser.add_argument('infile', help='Input file, for example \
    transposed_report.json')
parser.add_argument('outfile', help='Output file, chose a \
    file name, for example report.tsv')

args = parser.parse_args()

report = pd.read_json(args.infile)

report.to_csv(args.outfile, sep='\t', index=False)


print(f'Converted {args.infile} and saved copy to {args.outfile}')
