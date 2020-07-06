#!/usr/bin/env python

'''
    Usage:   python results_alleles_to_json.py <infile> <outfile>
'''

import pandas
import argparse

parser = argparse.ArgumentParser(description='This script converts a \
    tab-separated report file (.tsv) from chewBBACA \
    to json format.')

parser.add_argument('infile', help='Input file, for example \
    results_alleles.tsv')
parser.add_argument('outfile', help='Output file, chose a \
    file name, for example report_chewbbaca.json')

args = parser.parse_args()

report = pandas.read_csv(args.infile, sep='\t', index_col=0)
report.to_json(args.outfile, orient='index', indent=4)

print('Converted {} and saved copy to {}'.format(args.infile, args.outfile))
