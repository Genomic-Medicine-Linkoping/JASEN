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
aMRFinderPlus = $(aMRFinderPlus_find_pattern)

# ariba resistance find
# local "coding" db
ariba_resistance_out_dir = ariba_resistance
ariba_resistance_find_pattern_local_coding = motif_report_local.tsv
ariba_resistance_local_coding_out_file_end = motif_report_local_coding.tsv
# local "noncoding" db
ariba_resistance_find_pattern_local_noncoding = motif_report_nonc.tsv
ariba_resistance_local_noncoding_out_file_end = motif_report_local_noncoding.tsv
# joined local resistance hits
ariba_local_resistance = ariba_local_resistance.tsv
ariba_resfinder_resistance = ariba_resfinder_resistance.tsv


# resfinder db
ariba_resistance_find_pattern_resfinder = motif_report.tsv
ariba_resistance_resfinder_out_file_end = motif_report_local_resfinder.tsv


# mlst
mlst_out_dir = mlst
# ariba mlst
ariba_mlst_find_pattern = mlst_report.tsv
# ariba mlst joined output file
ariba_mlst = ariba_mlst.tsv

# tseemann mlst
tseemann_mlst_find_pattern = mlst.json
# tseemann mlst joined output file
tseemann_mlst = tseemann_mlst.tsv

# SpaTyper
spaTyper_out_dir = spaTyper
spaTyper = spaTyper.tsv

# Allowed number of threads
num_threads = 50

# Validation samples file
validation = validation_samples.tsv

# Final file with everything joined into it
summary = summary.tsv

define remove_ununique_adjacent_lines
	uniq $(1) > $(2)
	rm -f $(1)
endef

define remove_all_but_first_header
	head -n1 $(1) > temp.tsv
	sed '/^sample/d' $(1) >> temp.tsv
	$(call remove_ununique_adjacent_lines,temp.tsv,$(1))
endef


.PHONY: clean, gather_files, all

all: clean gather_files convert_json_to_tsv add_sample_species_columns subset_and_join_cols_of_interest join_to_validation_data

## gather_files:
## Note find "-L"
gather_files:
	./fetch_and_rename.sh $(results_file_root) $(aMRFinderPlus_find_pattern) 10 $(aMRFinderPlus_out_dir) $(aMRFinderPlus_find_pattern)
	./fetch_and_rename.sh $(results_file_root) $(ariba_resistance_find_pattern_local_coding) 10 $(ariba_resistance_out_dir) $(ariba_resistance_local_coding_out_file_end)
	./fetch_and_rename.sh $(results_file_root) $(ariba_resistance_find_pattern_local_noncoding) 10 $(ariba_resistance_out_dir) $(ariba_resistance_local_noncoding_out_file_end)
	./fetch_and_rename.sh $(results_file_root) $(ariba_resistance_find_pattern_resfinder) 10 $(ariba_resistance_out_dir) $(ariba_resistance_resfinder_out_file_end)
	./fetch_and_rename.sh $(results_file_root) $(ariba_mlst_find_pattern) 10 $(mlst_out_dir) $(ariba_mlst_find_pattern)
	./fetch_and_rename.sh $(results_file_root) $(tseemann_mlst_find_pattern) 10 $(mlst_out_dir) $(tseemann_mlst_find_pattern)
	mkdir -p $(spaTyper_out_dir) && \
	find -L $(results_file_root) -name "spa_*.tsv" -type f | parallel -j$(num_threads) "cp {} $(spaTyper_out_dir)/"
	find $(spaTyper_out_dir)/ -name "spa_*.tsv" -type f | parallel -j$(num_threads) "mv {} {= s:spa_:: =}"

## convert_json_to_tsv:
convert_json_to_tsv:
	$(CONDA_ACTIVATE) ; \
	find $(mlst_out_dir) -name "*.json" -type f | parallel -j$(num_threads) 'python json_to_tsv.py {} {.}.tsv'

## add_sample_species_columns:
add_sample_species_columns:
	$(CONDA_ACTIVATE) ; \
	find $(aMRFinderPlus_out_dir)/ $(ariba_resistance_out_dir)/ $(mlst_out_dir)/ $(spaTyper_out_dir)/ -name "*.tsv" -type f | parallel -j$(num_threads) './prepend_sample_species_cols.awk {} > {.}.prepend.tsv'

## subset_and_join_cols_of_interest: gather all results
subset_and_join_cols_of_interest:
	$(CONDA_ACTIVATE) ; \
	find $(aMRFinderPlus_out_dir)/ -name "*_aMRFinderPlus.prepend.tsv" -type f -print0 | xargs -0 -n1 -P$(num_threads) -I {} cut -d $$'\t' -f 1,2,8 {} > $(aMRFinderPlus)
	find $(mlst_out_dir)/ -name "*_mlst.prepend.tsv" -type f -print0 | xargs -0 -n1 -P$(num_threads) -I {} ./subset_ColsofInterest.py -m {} >> $(tseemann_mlst)
	find $(mlst_out_dir)/ -name "*_report.prepend.tsv" -type f -print0 | xargs -0 -I {} -P$(num_threads) cut -d $$'\t' -f 1-3 {} > $(ariba_mlst)
	find $(ariba_resistance_out_dir)/ -name "*_local*coding.prepend.tsv" -type f -print0 | xargs -0 -I {} -P$(num_threads) cut -d $$'\t' -f 1,2,4 {} > $(ariba_local_resistance)
	find $(ariba_resistance_out_dir)/ -name "*_motif_report_local_resfinder.prepend.tsv" -type f -print0 | xargs -0 -I {} -P$(num_threads) cut -d $$'\t' -f 1,2,4 {} > $(ariba_resfinder_resistance)
	find $(spaTyper_out_dir)/ -name "*.prepend.tsv" -type f -print0 | xargs -0 -I {} -P$(num_threads) cut -d $$'\t' -f 1,2,5 {} | sed 's/.tsv//' > $(spaTyper)
	$(call remove_all_but_first_header,$(aMRFinderPlus))
	$(call remove_all_but_first_header,$(tseemann_mlst))
	$(call remove_all_but_first_header,$(ariba_mlst))
	$(call remove_all_but_first_header,$(ariba_local_resistance))
	$(call remove_all_but_first_header,$(ariba_resfinder_resistance))
	$(call remove_all_but_first_header,$(spaTyper))

## join_to_validation_data: join all results to validation data
join_to_validation_data:
	$(CONDA_ACTIVATE) ; \
	jupyter nbconvert --inplace --to notebook --execute join_and_collapse.ipynb

## clean: Remove sub-directories and intermediary files
clean:
	rm -rf $(aMRFinderPlus_out_dir) ; \
	rm -rf $(ariba_resistance_out_dir) ; \
	rm -rf $(mlst_out_dir) ; \
	rm -rf $(spaTyper_out_dir) ; \
	rm -f $(tseemann_mlst) : \
	rm -f temp.tsv ; \
	rm -f all_data.tsv

## help: Show this message
help:
	@grep '^##' ./Makefile