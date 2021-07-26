SHELL = /bin/bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

.PHONY: all, \
run, \
clear_files, \
run_samples, \
update_subm, \
build_containers, \
preprocess_genomes

# Commands necessary for using conda env:s
CURRENT_CONDA_ENV_NAME = nf
# https://stackoverflow.com/a/55696820
# Note that the extra activate is needed to ensure that the activate floats env to the front of PATH
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate $(CURRENT_CONDA_ENV_NAME)

RG = assets/ref_genomes
PT = assets/prodigal_training_files

CONT_NAME = jasen_2021-07-07.sif
PROJECT_ROOT = /home/rada/Documents/CGL/JASEN

WORKDIR = $(PROJECT_ROOT)/work
IMAGE = $(PROJECT_ROOT)/container/$(CONT_NAME)

# Name of the species profile
SPECIES = Staphylococcus_aureus
# To which directory inside work/results should the output files come?
SAMPLE_ID = Staphylococcus_saprophyticus_Stam-121

# Command for running one sample
RUN = nextflow run main.nf -profile local,singularity,$(SPECIES) -resume --sample_ID $(SAMPLE_ID)

# Git branches
UPSTR_NAME = origin
UPSTR_BRANCH = main
CURR_BRANCH = local_aribadb


all: clear_files preprocess_genomes

clear_files:
	@echo ""
	@echo "Raderar både nedladdade och träningsfilerna:"
	@echo ""
	@(rm -f $(RG)/*.gz $(RG)/*.fna $(PT)/*.trn)
	@echo ""

build_containers:
	cd container; \
	sudo -E singularity build jasen_`date +%Y-%m-%d`.sif Singularity && \
	sudo -E /usr/local/bin/singularity build jasen_tidyverse_`date +%Y-%m-%d`.sif Singularity_tidyverse

preprocess_genomes:
	rm -f assets/ref_genomes/md5sums.txt
	bash preprocess_genomes.sh

run:
	@mkdir -p work
	rm -rf assets/references
	$(CONDA_ACTIVATE) ; \
	$(RUN)

update_subm:
	cd assets/var-genes-ro ; \
	# /usr/bin/git submodule update --remote --merge ; \
	/usr/bin/git fetch $(UPSTR_NAME) ; \
	/usr/bin/git merge $(UPSTR_NAME)/$(UPSTR_BRANCH) ; \
	cd .. ; \
	/usr/bin/git status ; \
	/usr/bin/git commit -am "Update submodule"
	# TODO: Could add the rest of the steps here as well
	# /usr/bin/git push $(UPSTR_NAME) $(CURR_BRANCH)

run_samples:
	$(CONDA_ACTIVATE) ; \
	bash run_samples.sh