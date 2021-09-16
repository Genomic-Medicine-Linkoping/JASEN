SHELL = /bin/bash
.ONESHELL:
.SHELLFLAGS := -e -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules


# Commands necessary for using conda env:s
CURRENT_CONDA_ENV_NAME = make-validation-summary
# https://stackoverflow.com/a/55696820
# Note that the extra activate is needed to ensure that the activate floats env to the front of PATH
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate $(CURRENT_CONDA_ENV_NAME)

results_file_root = results/

# aMRFinderPlus resistance find
aMRFinderPlus_out_dir = aMRFinderPlus
aMRFinderPlus_find_pattern = aMRFinderPlus.tsv

# ariba resistance find
# local "coding" db
ariba_resistance_out_dir = ariba_resistance
ariba_resistance_find_pattern_local_coding = motif_report_local.tsv
ariba_resistance_local_coding_out_file_end = motif_report_local_coding.tsv
# local "noncoding" db
ariba_resistance_find_pattern_local_noncoding = motif_report_nonc.tsv
ariba_resistance_local_noncoding_out_file_end = motif_report_local_noncoding.tsv
# resfinder db
ariba_resistance_find_pattern_resfinder = motif_report.tsv
ariba_resistance_resfinder_out_file_end = motif_report_local_resfinder.tsv

# mlst
mlst_out_dir = mlst
# ariba mlst
ariba_mlst_find_pattern = mlst_report.tsv
# tseemann mlst
tseemann_mlst_find_pattern = mlst.json

# SpaTyper
spaTyper_out_dir = spaTyper


.PHONY: clean, gather_files, all

all: clean gather_files 

# Note find "-L"
gather_files:
	./fetch_and_rename.sh $(results_file_root) $(aMRFinderPlus_find_pattern) 10 $(aMRFinderPlus_out_dir) $(aMRFinderPlus_find_pattern)
	./fetch_and_rename.sh $(results_file_root) $(ariba_resistance_find_pattern_local_coding) 10 $(ariba_resistance_out_dir) $(ariba_resistance_local_coding_out_file_end)
	./fetch_and_rename.sh $(results_file_root) $(ariba_resistance_find_pattern_local_noncoding) 10 $(ariba_resistance_out_dir) $(ariba_resistance_local_noncoding_out_file_end)
	./fetch_and_rename.sh $(results_file_root) $(ariba_resistance_find_pattern_resfinder) 10 $(ariba_resistance_out_dir) $(ariba_resistance_resfinder_out_file_end)
	./fetch_and_rename.sh $(results_file_root) $(ariba_mlst_find_pattern) 10 $(mlst_out_dir) $(ariba_mlst_find_pattern)
	./fetch_and_rename.sh $(results_file_root) $(tseemann_mlst_find_pattern) 10 $(mlst_out_dir) $(tseemann_mlst_find_pattern)
	mkdir -p $(spaTyper_out_dir) && \
	find -L $(results_file_root) -name "spa_*.tsv" -type f | parallel -j+0 "cp {} $(spaTyper_out_dir)/"
	find $(spaTyper_out_dir)/ -name "spa_*.tsv" -type f | parallel -j+0 "mv {} {= s:spa_:: =}"

convert_json_to_tsv:
	$(CONDA_ACTIVATE) ; \
	find $(mlst_out_dir) -name "*.json" -type f | parallel -j+0 'python json_to_tsv.py {} {.}.tsv'

add_sample_species_columns:
	$(CONDA_ACTIVATE) ; \
	find $(aMRFinderPlus_out_dir)/ $(ariba_resistance_out_dir)/ $(mlst_out_dir)/ $(spaTyper_out_dir)/ -name "*.tsv" -type f | parallel -j+0 './prepend_sample_species_cols.awk {} > {.}.prepend.tsv'

clean:
	rm -rf $(aMRFinderPlus_out_dir) ; \
	rm -rf $(ariba_resistance_out_dir) ; \
	rm -rf $(mlst_out_dir) ; \
	rm -rf $(spaTyper_out_dir)