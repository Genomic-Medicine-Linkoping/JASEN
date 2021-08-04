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
	spa_dir = Path(args.spa_dir)
	ariba_mlst_dir = Path(args.ariba_mlst_dir)
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
		genes = set(df["ref_name"].to_list())
	if (any(nt.species_id==species_id for nt in species_id_res_genes_list)):
		# Obtain already in list stored NamedTuple instance
		current_nt = [nt for nt in species_id_res_genes_list if nt.species_id==species_id][0]
		# Use sets to store gene names so that only unique genes are stored
		current_nt.genes.update(genes)
	else:
		# If the NamedTuple with specific species_id doesn't exist in the list yet,
		# add it to the list with genes inside a set
			species_id_res_genes_list.append(Resistance_genes(species_id, genes))

	spa_types = {}

	for path in sorted(spa_dir.glob(args.spa_typer_pattern)):
		df = pd.read_csv(path, sep='\t')
		# Extract species_id data from file name, e.g. Staphylococcus_aureus_SASPA17S-03
		species_id = re.sub(r'spa_(.+)', r'\1', str(path.stem)) 
		# Extract spa-type, e.g. t1339
		spa_type = df.iloc[0,2]
		spa_types[species_id] = spa_type

	ariba_mlst = {}

	for path in sorted(ariba_mlst_dir.glob(args.ariba_mlst_pattern)):
		df = pd.read_csv(path, sep='\t')
		# Extract species_id data from file name, e.g. Staphylococcus_aureus_SASPA17S-03
		species_id = re.sub('(.+)_mlst_report', r'\1', str(path.stem)) 
		# Extract ariba mlst-type, e.g. t1339
		ariba_mlst_type = df.iloc[0,0]
		ariba_mlst[species_id] = ariba_mlst_type

with open(tsv, 'w') as f:
	# Write tsv header
	f.write(f'Prov-nummer\tArt\tKända resistensgener samt PVL\tKänd MLST eller Spa-typ\n')
	for current_si in species_id_res_genes_list:
		mlst_type = st.get(current_si.species_id)
		si = current_si.species_id.split('_')
		species = si[0][0] + '. ' + si[1]
		id = si[2:][0]
			# Extract spa-type if it exists
			spa_type = spa_types.get(current_si_nt.species_id)
			# Extract ariba mlst-type
			ariba_mlst_type = ariba_mlst.get(current_si_nt.species_id)
			
		# Write nothing if MLST type is missing
		if (mlst_type=='-'):
				# Do we write spa-type or not when MLST type is missing
				if (not spa_type):
					f.write(f'{id}\t{species}\t{genes_str}\taST{ariba_mlst_type}\n')
		else:
					f.write(f'{id}\t{species}\t{genes_str}\taST{ariba_mlst_type}, {spa_type}\n')
			else:
				# Do we write spa-type or not when MLST type is found
				if (not spa_type):
					f.write(f'{id}\t{species}\t{genes_str}\taST{ariba_mlst_type}, ST{mlst_type}\n')
				else:
					f.write(f'{id}\t{species}\t{genes_str}\taST{ariba_mlst_type}, ST{mlst_type}, {spa_type}\n')



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
