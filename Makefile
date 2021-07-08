# https://stackoverflow.com/a/55696820
SHELL = /bin/bash

CURRENT_CONDA_ENV_NAME = nf
# Note that the extra activate is needed to ensure that the activate floats env to the front of PATH
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate $(CURRENT_CONDA_ENV_NAME)


#FASTQS?=.github/data/fastqs/
RG = assets/ref_genomes
PT = assets/prodigal_training_files

CONT_NAME = jasen_2021-07-07.sif
PROJECT_ROOT = /home/rada/Documents/CGL/JASEN

WORKDIR = $(PROJECT_ROOT)/work
IMAGE = $(PROJECT_ROOT)/container/$(CONT_NAME)

# Name of the species profile
SPECIES = Staphylococcus_aureus
# SPECIES = Escherichia_coli
# To which directory inside work/results should the output files come?
# SAMPLE_ID = Klebsiella_pneumoniae_p1
# SAMPLE_ID = Escherichia_coli_p1
#SAMPLE_ID = Staphylococcus_aureus_prov1
SAMPLE_ID = Staphylococcus_aureus_prov2

# env TZ="Europe/Stockholm" and -B /run 
# fixes a problem described in: https://github.com/truatpasteurdotfr/singularity-docker-fedora30-brave/issues/3
#RUN = env TZ="Europe/Stockholm" /usr/local/bin/singularity exec -B /run -B $(PROJECT_ROOT):/external -B $(WORKDIR):/out $(IMAGE) nextflow -C /external/nextflow.config run main.nf -profile local,singularity,$(SPECIES) -resume
RUN = nextflow run main.nf -profile local,singularity,$(SPECIES) -resume --sample_ID $(SAMPLE_ID)

UPSTR_NAME = origin
UPSTR_BRANCH = main
CURR_BRANCH = local_aribadb

all: clear_files download_bacterial_genomes create_prodigal_trn_files uncompress_genomes

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
	bash download_genomes.sh

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
