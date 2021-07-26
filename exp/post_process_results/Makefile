SHELL = /bin/bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

# Scripts
COMP_MLST = compare_mlst.py
COMP_RESIST = compare_resistence.py
GATHER = gather_files.py
JOIN = join_nf_results.py

# Gathered results dirs
MLST = /home/rada/Documents/CGL/JASEN/exp/post_process_results/mlst_json
RESISTANCE = /home/rada/Documents/CGL/JASEN/exp/post_process_results/resistens

# Where to look for results
INDIR = /data/CGL/JASEN/results

VAL = validerings_prov.tsv
JOINED_RES = joined.tsv

RES_OUT = resistance_merged.tsv
MLST_OUT = mlst_merged.tsv

# Dummies
GATHERED = gather.rdy

.PHONY: all, clean

all: $(RES_OUT) $(MLST_OUT) #$(GATHERED) $(JOINED_RES)

# Create a dummy file to keep track of having gathered the results of interest into the right directories
$(GATHERED): 
	python3 $(GATHER) $(INDIR) $(MLST)
	python3 $(GATHER) $(INDIR) $(RESISTANCE)	
	touch $@

# Join all mlst and resistance JASEN results into one file
$(JOINED_RES): $(JOIN) $(GATHERED)
	python3 $< $(MLST) $(RESISTANCE) $@

# Create final output files: comparison of validation and JASEN resistance findings
$(RES_OUT): $(COMP_RESIST) $(VAL) $(JOINED_RES)
	python3 $^ $@

# Create final output files: comparison of validation and JASEN MLST findings
$(MLST_OUT): $(COMP_MLST) $(VAL) $(JOINED_RES)
	python3 $^ $@

clean:
	rm -rf $(RES_OUT) $(MLST_OUT) $(GATHERED) $(JOINED_RES)