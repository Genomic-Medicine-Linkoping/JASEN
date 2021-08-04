SHELL := /bin/bash
.ONESHELL:
#.SHELLFLAGS := -eu -o pipefail -c
.SHELLFLAGS := -e -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

.PHONY: all, \
run, \
clean, \
run_samples, \
update_subm, \
preprocess

# Commands necessary for using conda env:s
CURRENT_CONDA_ENV_NAME = nf
# https://stackoverflow.com/a/55696820
# Note that the extra activate is needed to ensure that the activate floats env to the front of PATH
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate $(CURRENT_CONDA_ENV_NAME)

RG = assets/ref_genomes
PT = assets/prodigal_training_files

CONT_NAME = library://ljmesi/jasen/main.sif
PROJECT_ROOT = $(PWD)

WORKDIR = $(PROJECT_ROOT)/work
IMAGE = $(PROJECT_ROOT)/container/$(CONT_NAME)

# Name of the species profile
SPECIES = Staphylococcus_saprophyticus
# This is the name of the input directory (with the fastq files) in assets/sequencing_data
# and also the name of the output directory in results/
SAMPLE_ID = Staphylococcus_saprophyticus_Stam-121

# Command for running one sample
RUN = nextflow run main.nf -profile local,singularity,$(SPECIES) -resume --sample_ID $(SAMPLE_ID) $(ARGS)

# Git branches
UPSTR_NAME = origin
UPSTR_BRANCH = main
CURR_BRANCH = ro-implementation

## all: Run by default one sample
all: run

## clean: Remove all downloaded genome files, prodigal training files and checksum file
clean:
	@echo ""
	@echo "Remove all downloaded genome files, prodigal training files and checksum file"
	@echo ""
	rm -f $(RG)/*.gz $(RG)/*.fna $(PT)/*.trn && \
	rm -f $(RG)/md5sums.txt
	@echo ""

## preprocess: Download, uncompress and create prodigal training files of genomes and create md5sum:s
preprocess: clean
	cp bin/preprocess_genomes.sh assets/genome_data.tsv . && \
	bash preprocess_genomes.sh $(CONT_NAME) && \
	rm -f preprocess_genomes.sh genome_data.tsv

## run: Run one sample with $(SPECIES) and $(SAMPLE_ID) located in assets/sequencing_data
run:
	@mkdir -p work
	rm -rf assets/references
	$(CONDA_ACTIVATE) ; \
	$(RUN)

## update_subm: Update assets/var-genes-ro submodule and push it to ro-implementation remote branch
update_subm:
	cd assets/var-genes-ro ; \
	# /usr/bin/git submodule update --remote --merge ; \
	/usr/bin/git fetch $(UPSTR_NAME) ; \
	/usr/bin/git merge $(UPSTR_NAME)/$(UPSTR_BRANCH) ; \
	cd .. ; \
	/usr/bin/git status ; \
	/usr/bin/git add var-genes-ro ; \
	/usr/bin/git commit -m "Update submodule" ; \
	/usr/bin/git push $(UPSTR_NAME) $(CURR_BRANCH)

## run_samples: Run pipeline with a set of samples in: assets/sequencing_data
run_samples:
	cp bin/run_samples.sh . && \
	$(CONDA_ACTIVATE) ; \
	bash run_samples.sh && \
	rm -f run_samples.sh

## help: show this message.
help:
	@grep '^##' ./Makefile