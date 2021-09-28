SHELL := /bin/bash
.ONESHELL:
#.SHELLFLAGS := -eu -o pipefail -c
.SHELLFLAGS := -e -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

.PHONY: all, \
run, \
run_tax_analysis, \
clean, \
run_samples, \
update_subm, \
preprocess, \
download_latest_spa_db, \
help

# Commands necessary for using conda env:s
CURRENT_CONDA_ENV_NAME = nf
# https://stackoverflow.com/a/55696820
# Note that the extra activate is needed to ensure that the activate floats env to the front of PATH
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate $(CURRENT_CONDA_ENV_NAME)

RG = assets/ref_genomes
PT = assets/prodigal_training_files

CONT_NAME = library://ljmesi/jasen/main.sif
PROJECT_ROOT = $(PWD)

IMAGE = $(PROJECT_ROOT)/container/$(CONT_NAME)

SPA_DB = assets/spa-typing
KAIJU_DB = /data/CGL/JASEN/kaiju-db
KRAKEN_DB_DIR = /data/CGL/JASEN
KRAKEN_DB_NAME = K2DB

# Name of the species profile
SPECIES = Staphylococcus_aureus
# This is the name of the input directory (with the fastq files) in assets/sequencing_data
# and also the name of the output directory in results/

RESULTS_DIR = $(PROJECT_ROOT)/results/$(SAMPLE_ID)

# Command for running one sample
RUN = nextflow run main.nf -profile local,singularity,$(SPECIES) -resume --sample_ID $(SAMPLE_ID) -with-report $(RESULTS_DIR)/report.html -with-timeline $(RESULTS_DIR)/timeline.html -with-dag $(RESULTS_DIR)/dag.html $(ARGS)

# Git branches
UPSTR_NAME = origin
UPSTR_BRANCH = main
CURR_BRANCH = ro-implementation

## all: Run by default one sample
all: run

## PREPARING/UPDATING THE PIPELINE
## preprocess: Download, uncompress and create prodigal training files of genomes and create md5sum:s
preprocess: clean download_latest_spa_db
	cp bin/preprocess_genomes.sh assets/genome_data.tsv . && \
	bash preprocess_genomes.sh $(CONT_NAME) && \
	rm -f preprocess_genomes.sh genome_data.tsv

## create_kaijudb: Create database that kaiju uses in its taxonomic assignments
create_kaijudb:
	$(CONDA_ACTIVATE) ; \
	mkdir -p $(KAIJU_DB) && \
	cd $(KAIJU_DB) && \
	kaiju-makedb -t 50 -s nr_euk

## create_kraken2db: Create kraken2 standard DB with: NCBI taxonomic information, the complete genomes in RefSeq for 
## the bacterial, archaeal, and viral domains, along with the human genome and a collection of known vectors (UniVec_Core)
##
create_kraken2db:
	$(CONDA_ACTIVATE) ; \
	mkdir -p $(KRAKEN_DB_DIR) && \
	cd $(KRAKEN_DB_DIR) && \
	kraken2-build --standard --threads 50 --use-ftp --db $(KRAKEN_DB_NAME)

## download_latest_spa_db: Download the newest sparepeats.fasta and spatypes.txt files to assets/spa-typing
download_latest_spa_db:
	curl -o $(SPA_DB)/sparepeats.fasta https://spa.ridom.de/dynamic/sparepeats.fasta && \
	curl -o $(SPA_DB)/spatypes.txt https://spa.ridom.de/dynamic/spatypes.txt

## RUNNING THE PIPELINE
## run_tax_analysis: Run a main pipeline preceding taxonomic analysis of the sample fastq files
run_tax_analysis:
	rm -rf work
	@mkdir -p work
	rm -rf assets/references
	$(CONDA_ACTIVATE) ; \
	$(RUN) --run_tax_analysis true $(ARGS)

## run: Run one sample with $(SPECIES) and $(SAMPLE_ID) located in assets/sequencing_data
## NB: Remember to run make run_tax_analysis before running this
##
run:
	@mkdir -p work
	$(CONDA_ACTIVATE) ; \
	$(RUN)

## run_samples: Run pipeline with a set of samples in: assets/sequencing_data
run_samples:
	cp bin/run_samples.sh . && \
	$(CONDA_ACTIVATE) ; \
	bash run_samples.sh && \
	rm -f run_samples.sh

## UTILITIES
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

## build_main_singularity: (Re)build singularity image
build_main_singularity:
	cd container ; \
	rm -f main.sif && \
	sudo singularity build main.sif Singularity

## push_to_cloud: Sign and push built image to Sylabs cloud
## NB: Remember to rename on the cloud existing image to something else than 'latest' before running this
##
push_to_cloud:
	cd container ; \
	singularity sign main.sif && \
	singularity push main.sif library://ljmesi/jasen/main.sif:latest

## clean: Remove all downloaded genome files, prodigal training files and checksum file
clean:
	@echo ""
	@echo "Remove all downloaded genome files, prodigal training files and checksum file"
	@echo ""
	rm -f $(RG)/*.gz $(RG)/*.fna $(PT)/*.trn && \
	rm -f $(RG)/md5sums.txt
	@echo ""

## help: Show this message
help:
	@grep '^##' ./Makefile
