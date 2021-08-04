#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""This script joins separate tab separated files in various directories into one tsv file."""

import argparse
from pathlib import Path
import pandas as pd
import re
from collections import namedtuple

def main(args):
mlst_dir = Path(args.mlst_dir)
resistance_dir = Path(args.resistance_dir)
tsv = Path(args.outfile)
json_pat = args.json_pattern

st = {}

# Store all found MLST types in a dict 
for path in sorted(mlst_dir.glob(json_pat)):
	df = pd.read_json(path)
	sequence_type = df.iloc[0]['sequence_type']
	id = re.sub(json_pat[1:], '', str(path.name))
	st[id] = str(sequence_type)

# Define the namedtuple
Resistance_genes = namedtuple('Resistance_genes', ['species_id','genes'])
# List for storing species_id specific NamedTuples  
species_id_res_genes_list = []

# 
for path in sorted(resistance_dir.glob(args.res_pattern)):
	df = pd.read_csv(path, sep='\t')
	# Extract species and ID data from file name
	species_id = re.sub(r'_motif_report.*\.tsv', '', str(path.name))
	# Convert all genes in the current dataframe to a list
	genes = df["ref_name"].to_list()
	# Check if the list contains already the current species_id inside a namedtuple
	if (any(nt.species_id==species_id for nt in species_id_res_genes_list)):
		# Obtain already in list stored NamedTuple instance
		current_nt = [nt for nt in species_id_res_genes_list if nt.species_id==species_id][0]
		# Use sets to store gene names so that only unique genes are stored
		current_nt.genes.update(genes)
	else:
		# If the NamedTuple with specific species_id doesn't exist in the list yet,
		# add it to the list with genes inside a set
		species_id_res_genes_list.append(Resistance_genes(species_id, set(genes)))
	# break



with open(tsv, 'w') as f:
	# Write tsv header
	f.write(f'Prov-nummer\tArt\tKända resistensgener samt PVL\tKänd MLST eller Spa-typ\n')
	for current_si in species_id_res_genes_list:
		mlst_type = st.get(current_si.species_id)
		si = current_si.species_id.split('_')
		species = si[0][0] + '. ' + si[1]
		id = si[2:][0]
		genes_str = ', '.join(current_si.genes)
		# Write nothing if MLST type is missing
		if (mlst_type=='-'):
			f.write(f'{id}\t{species}\t{genes_str}\t\n')
		else:


if __name__ == '__main__':
	# Create the parser
	parser = argparse.ArgumentParser(prog='join_nf_results.py',
									 allow_abbrev=False,
									 usage='%(prog)s [options] <mlst_dir> <resistance_dir> <outfile>',
									 description='',
									 epilog='---')

	parser.add_argument('mlst_dir',
			type=str,
			help='Place where all the mlst typing results are')

	parser.add_argument('resistance_dir',
			type=str,
			help='Place where all the ariba resistance gene search results are')

	parser.add_argument('spa_dir',
			type=str,
			help='Place where all the spaTyper results are')
	
	parser.add_argument('ariba_mlst_dir',
			type=str,
			help='Place where all the ariba mlst results are')

	parser.add_argument('outfile', 
			type=str,
			help='Place where we want our joined jsons be saved as one tsv file')

	parser.add_argument('--json-pattern',
			'-j',
			type=str,
			default='*_mlst.json',
			help='Globbing pattern to search mlst json files with')

	parser.add_argument('--res-pattern',
			'-p',
			type=str,
			default='*.tsv',
			help='Globbing pattern to search ariba tsv files with')

	parser.add_argument('--spa-typer-pattern',
			'-s',
			type=str,
			default='spa_*.tsv',
			help='Globbing pattern to search spaTyper tsv files with')

	parser.add_argument('--ariba-mlst-pattern',
			'-a',
			type=str,
			default='*_mlst_report.tsv',
			help='Globbing pattern to search ariba mlst tsv files with')

	# Execute the parse_args() method
	args = parser.parse_args()
	main(args)
