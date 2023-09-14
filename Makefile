# ==============================================================================
# Brief explanation of Makefile syntax
# ==============================================================================
# Since Makefiles are a bit special, and there are not a lot of easy to read
# tutorials out there, here follows a very brief introduction to the syntax in
# Makefiles.
#
# Basics
# ------------------------------------------------------------------------------
#
# Makefiles are in many ways similar to bash-scripts, but unlike bash-scripts
# they are not written in a linear sequential fashion, but rather divided into
# so called rules, that are typically tightly connected to specific output file
# paths or file path pattern.
#
# Firstly, Makefiles should be named `Makefile` and be put inside a folder
# where one wants to run the make command.
#
# The rule syntax
# ------------------------------------------------------------------------------
#
# The syntax of a rule is at its core very simple:
#
# output file(s)/target(s) : any dependent input files
#     commands to produce the output files from input files
#     possibly more commands
#
# So, the significant part of the rule syntax is the ':'-character, as well as
# the indentation of the rows following the head of the rule, to indicate the
# "recipe" or the commands to produce the outputs of the rule.
#
# A rule can also have just a name, and no concrete output file. That is, it
# would have the form:
#
# name_of_the_rule : any dependent input files
#     one commands
#     more commands
#
# Now, there is one big caveat here, related to our scripts: Make will rebuild
# a target as soon as any of its inputs are updated, or have a newer timestamp
# than the target. This is typically not desired for us, since we might have
# files unpacked with zip with all kinds of different dates, and for a one-run
# installation, we are mostly interested in wheter an output file already
# exists or not, not so much about timestamps.
#
# To change so that Make only cares about whether files exist, and not timestamps,
# one can add a | character before those input files, like so:
#
# name_of_the_rule : input files where timestamp matters | input files where only existence matters
#     one commands
#     more commands
#
# Of course, one can have everything on the right side of the | character, so:
#
# name_of_the_rule : | input files where only existence matters
#     one commands
#     more commands
#
# Running Makefiles
# ------------------------------------------------------------------------------
#
# To run a named rule, you would then do:
#
# $ make name_of_the_rule
#
# For normal rules, you would hit:
#
# $ make <outputfile>
#
# Tip: Type "make" and hit tab twice in the shell, to show available targets to
# run in the current makefile.
#
# Special variables and more
# ------------------------------------------------------------------------------
#
# Inside the commands of the rules, one can write pretty much bash, especially
# if setting `SHELL := /bin/bash` in the beginning of the Makefile which we have
# done below.
#
# There are a some differences though:
#
# 1. Variable syntax using only a single $-character refer to MAKE-variables.
#    Thus, to use variables set in bash, you have to use double-$.
#    So:
#    echo $(this-is-a-make-variable) and $$(this-is-a-bash-variable)
#
#    This also goes for command substitution, so rather than:
#
#    	echo "Lines in file:" $(wc -l somefile.txt)
#
#    ... you would write:
#
#    	echo "Lines in file:" $$(wc -l somefile.txt)
#
# 2. Makefiles also use some special variables with special meaning, to access
#    things like the output and input files:
#
#    $@   The output target (In case of multiple targets, the same command will
#         be run once for every output file, feeding this variable with only a
#         single one at a time)
#
#    $<   The first input file or dependency
#
#    $^   ALL input files / dependencies.
#
# 3. Another speciality is that adding a '@' character before any command, will
#    stop this command from being printed out when it is executed (only its output
#    will be printed).
#
# That's all for now. For the official docs, see:
# https://www.gnu.org/software/make/manual/make.html
#
# Some more caveats
# ------------------------------------------------------------------------------
#  Here are some things that might not be obvious at first:
#
#  To create a rule with a single output but many dependencies, you can add
#  these dependencies on multiple lines by using the continuation character \
#  like so:
#
#  name_of_rule: dependency1 \
#  		dependency2 \
#  		dependency3 \
#  		dependency4 \

# ==============================================================================
# Various definitions
# ==============================================================================
# Make sure bash is used as shell, for consistency and to make some more
# advanced scripting possible than with /bin/sh
SHELL := /bin/bash

# Define path variables
SCRIPT_DIR := $(shell pwd)
ASSETS_DIR := $(shell realpath $(SCRIPT_DIR)/assets/)
CONTAINER_DIR := $(realpath $(SCRIPT_DIR)/container/)
PRODIGAL_TRAINING_DIR := $(ASSETS_DIR)/prodigal_training_files
# The root folder where the pipeline is currently located. To be mounted into
# the Singularity containers below.
MNT_ROOT := /$(shell readlink -f . | cut -d"/" -f2)
INSTALL_LOG := "$(SCRIPT_DIR)/.install.log"

# Declare that these commands don't produce outputs
.PHONY: print_vars install update_amrfinderplus update_mlst_db update_blast_db update_virulencefinder_db

define log_message
	@echo "--------------------------------------------------------------------------------" | tee -a $(INSTALL_LOG);
	@echo "$$(date "+%Y-%m-%d %H:%M:%S"): $1" | tee -a $(INSTALL_LOG);
	@echo "--------------------------------------------------------------------------------" | tee -a $(INSTALL_LOG);
endef

# ==============================================================================
# Main rules
# ==============================================================================

print_paths:
	@echo "SCRIPT_DIR:" $(SCRIPT_DIR)
	@echo "ASSETS_DIR:" $(ASSETS_DIR)
	@echo "CONTAINER_DIR:" $(CONTAINER_DIR)
	@echo "MNT_ROOT:" $(MNT_ROOT)

# NOTE: We don't run *_all for kpneumoniae, since prodigal training file is missing
install: build_containers \
	update_amrfinderplus \
	update_mlst_db \
	update_blast_db \
	update_virulencefinder_db \
	saureus_all \
	ecoli_all \
	kpneumoniae_unpack_cgmlst_schema \
	mtuberculosis_all

check:	check_chewbbaca \
	check_bwa \
	check_blastdb


# ==============================================================================
# Build containers
# ==============================================================================

build_containers:
	cd $(CONTAINER_DIR) \
	&& make all; \
	cd -

# ==============================================================================
# Update databases
# ==============================================================================

# -----------------------------
# Update AMRFinderPlus database
# -----------------------------
# For more info, see:
# https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/
AMRFINDERDB_DIR := $(ASSETS_DIR)/amrfinder_db

update_amrfinderplus: $(AMRFINDERDB_DIR)/latest

$(AMRFINDERDB_DIR)/latest:
	$(call log_message,"Starting update of AMRFinderPlus database ...")
	singularity exec \
		--bind $(MNT_ROOT) \
		$(CONTAINER_DIR)/amrfinderplus.sif \
		amrfinder_update \
		--database $(AMRFINDERDB_DIR) |& tee -a $(INSTALL_LOG)

# -----------------------------
# Update MLST database
# -----------------------------
MLSTDB_DIR := $(ASSETS_DIR)/mlst_db

update_mlst_db: $(MLSTDB_DIR)/pubmlst/dbases.xml

$(ASSETS_DIR)/mlst_db/pubmlst/dbases.xml:
	$(call log_message,"Starting update of MLST database ...")
	cd $(MLSTDB_DIR) \
	&& bash mlst-download_pub_mlst.sh |& tee -a $(INSTALL_LOG)

# -----------------------------
# Update Blast database
# -----------------------------
update_blast_db: $(MLSTDB_DIR)/blast

$(MLSTDB_DIR)/blast:
	$(call log_message,"Starting update of Blast database")
	cd $(MLSTDB_DIR) \
	&& singularity exec \
		--bind $(MNT_ROOT) \
		$(CONTAINER_DIR)/blast.sif \
		bash $(MLSTDB_DIR)/mlst-make_blast_db.sh |& tee -a $(INSTALL_LOG)

# -----------------------------
# Update VirulenceFinder database
# -----------------------------
VIRULENCEFINDERDB_DIR := $(ASSETS_DIR)/virulencefinder_db
KMA_DIR := $(ASSETS_DIR)/kma

update_virulencefinder_db: $(VIRULENCEFINDERDB_DIR)/stx.name

$(ASSETS_DIR)/virulencefinder_db/stx.name:
	$(call log_message,"Starting update of VirulenceFinder database")
	cd $(ASSETS_DIR)/kma \
	&& make \
	&& cd $(VIRULENCEFINDERDB_DIR) \
	&& export PATH=$(ASSETS_DIR)/kma:$$PATH \
	&& singularity exec	--bind $(MNT_ROOT) $(CONTAINER_DIR)/pythonScripts.sif \
		python3 INSTALL.py \
		$(KMA_DIR)/kma_index |& tee -a $(INSTALL_LOG) \
	&& cd $(ASSETS_DIR)/resfinder_db \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/pythonScripts.sif \
		python3 INSTALL.py \
		$(KMA_DIR)/kma_index |& tee -a $(INSTALL_LOG) \
	&& cd $(ASSETS_DIR)/pointfinder_db \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/pythonScripts.sif \
		python3 INSTALL.py \
		$(KMA_DIR)/kma_index |& tee -a $(INSTALL_LOG)

# ==============================================================================
# Download, index and prep reference genomes for organisms
# ==============================================================================

# -----------------------------
# S. Aureus
# -----------------------------
.PHONY: saureus_all \
	saureus_download_reference \
	saureus_index_reference \
	saureus_download_cgmlst_schema \
	saureus_unpack_cgmlst_schema \
	saureus_prep_cgmlst_schema

saureus_all: saureus_download_reference \
	saureus_index_reference \
	saureus_download_cgmlst_schema \
	saureus_unpack_cgmlst_schema \
	saureus_prep_cgmlst_schema

SAUR_GENOMES_DIR := $(ASSETS_DIR)/genomes/staphylococcus_aureus
SAUR_CGMLST_DIR := $(ASSETS_DIR)/cgmlst/staphylococcus_aureus
SAUR_REFSEQ_ACC := NC_002951.2


saureus_download_reference: $(SAUR_GENOMES_DIR)/$(SAUR_REFSEQ_ACC).fasta

$(SAUR_GENOMES_DIR)/$(SAUR_REFSEQ_ACC).fasta:
	$(call log_message,"Downloading S. Aureus reference genome ...")
	mkdir -p $(SAUR_GENOMES_DIR) \
	&& cd $(SCRIPT_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/pythonScripts.sif \
		python3 bin/download_ncbi.py \
		-i $(SAUR_REFSEQ_ACC) \
		-o $(SAUR_GENOMES_DIR) |& tee -a $(INSTALL_LOG)


saureus_index_reference: $(SAUR_GENOMES_DIR)/$(SAUR_REFSEQ_ACC).fasta.bwt

$(SAUR_GENOMES_DIR)/$(SAUR_REFSEQ_ACC).fasta.bwt: $(SAUR_GENOMES_DIR)/$(SAUR_REFSEQ_ACC).fasta
	$(call log_message,"Indexing S. Aureus reference genome ...")
	cd $(SAUR_GENOMES_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/bwakit.sif \
		bwa index $$(basename $<) |& tee -a $(INSTALL_LOG)


saureus_download_cgmlst_schema: $(SAUR_CGMLST_DIR)/alleles/cgmlst_141106.zip

$(SAUR_CGMLST_DIR)/alleles/cgmlst_141106.zip:
	$(call log_message,"Downloading S. Aureus cgMLST schema ...")
	mkdir -p $(SAUR_CGMLST_DIR)/alleles &> /dev/null \
	&& cd $(SAUR_CGMLST_DIR)/alleles \
	&& wget https://www.cgmlst.org/ncs/schema/141106/alleles/ \
		-O $@ \
		--no-check-certificate |& tee -a $(INSTALL_LOG)


saureus_unpack_cgmlst_schema: $(SAUR_CGMLST_DIR)/alleles/SACOL2694.fasta

$(SAUR_CGMLST_DIR)/alleles/SACOL2694.fasta: $(SAUR_CGMLST_DIR)/alleles/cgmlst_141106.zip
	$(call log_message,"Unpacking S. Aureus cgMLST schema ...")
	cd $$(dirname $<) \
		&& unzip -Duq $$(basename $<) |& tee -a $(INSTALL_LOG)


saureus_prep_cgmlst_schema: $(SAUR_CGMLST_DIR)/alleles_rereffed/Staphylococcus_aureus.trn

$(SAUR_CGMLST_DIR)/alleles_rereffed/Staphylococcus_aureus.trn: $(SAUR_CGMLST_DIR)/alleles/SACOL2694.fasta
	$(call log_message,"Prepping S. Aureus cgMLST schema ...")
	cd $(SAUR_CGMLST_DIR) \
	&& echo "WARNING! Prepping cgMLST schema. This takes a looong time. Put on some coffee" \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/chewbbaca.sif \
		chewie PrepExternalSchema \
		-i $(SAUR_CGMLST_DIR)/alleles \
		-o $(SAUR_CGMLST_DIR)/alleles_rereffed \
		--cpu 1 \
		--ptf $(PRODIGAL_TRAINING_DIR)/Staphylococcus_aureus.trn |& tee -a $(INSTALL_LOG)

# -----------------------------
# E. coli
# -----------------------------

.PHONY: ecoli_all \
	ecoli_download_reference \
	ecoli_index_reference \
	ecoli_download_wgmlst_schema \
	ecoli_download_cgmlst_schema \
	ecoli_prep_ecoli_cgmlst_schema

ecoli_all: ecoli_index_reference \
	ecoli_download_wgmlst_schema \
	ecoli_prep_ecoli_cgmlst_schema

ECOLI_GENOMES_DIR := $(ASSETS_DIR)/genomes/escherichia_coli
ECOLI_WGMLST_DIR := $(ASSETS_DIR)/wgmlst/escherichia_coli
ECOLI_CGMLST_DIR := $(ASSETS_DIR)/cgmlst/escherichia_coli
ECOLI_REFSEQ_ACC := NC_000913.3


ecoli_download_reference: $(ECOLI_GENOMES_DIR)/$(ECOLI_REFSEQ_ACC).fasta

$(ECOLI_GENOMES_DIR)/$(ECOLI_REFSEQ_ACC).fasta:
	$(call log_message,"Downloading E. coli genome ...")
	cd $(SCRIPT_DIR) \
	&& mkdir -p $(ECOLI_GENOMES_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/pythonScripts.sif \
		python3 bin/download_ncbi.py \
		-i $(ECOLI_REFSEQ_ACC) \
		-o $(ECOLI_GENOMES_DIR) |& tee -a $(INSTALL_LOG)


ecoli_index_reference: $(ECOLI_GENOMES_DIR)/$(ECOLI_REFSEQ_ACC).fasta.bwt

$(ECOLI_GENOMES_DIR)/$(ECOLI_REFSEQ_ACC).fasta.bwt: $(ECOLI_GENOMES_DIR)/$(ECOLI_REFSEQ_ACC).fasta
	$(call log_message,"Indexing E. coli genome ...")
	cd $(ECOLI_GENOMES_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/bwakit.sif \
		bwa index $$(basename $<) |& tee -a $(INSTALL_LOG)


# Download Ecoli wgmlst INNUENDO schema
ecoli_download_wgmlst_schema: $(ECOLI_WGMLST_DIR)/alleles/ecoli_INNUENDO_wgMLST/Escherichia_coli.trn

$(ECOLI_WGMLST_DIR)/alleles/ecoli_INNUENDO_wgMLST/Escherichia_coli.trn:
	$(call log_message,"Downloading E. coli wgMLST schema ...")
	mkdir -p $(ECOLI_WGMLST_DIR)/alleles &> /dev/null \
	&& cd $(ECOLI_WGMLST_DIR)/alleles \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/chewbbaca.sif \
		chewie DownloadSchema \
		-sp 5 \
		-sc 1 \
		-o $(ECOLI_WGMLST_DIR)/alleles \
		--latest |& tee -a $(INSTALL_LOG)


# Download Ecoli cgmlst cgmlst.org schema
ecoli_download_cgmlst_schema: $(ECOLI_CGMLST_DIR)/alleles/ecoli_cgmlst_alleles_5064703.zip

$(ECOLI_CGMLST_DIR)/alleles/ecoli_cgmlst_alleles_5064703.zip:
	$(call log_message,"Downloading E. coli cgMLST schema ...")
	mkdir -p $(ECOLI_CGMLST_DIR)/alleles &> /dev/null \
	&& cd $(ECOLI_CGMLST_DIR)/alleles \
	&& wget -O $$(basename $@) https://www.cgmlst.org/ncs/schema/5064703/alleles/ --no-check-certificate |& tee -a $(INSTALL_LOG)


# Unpack Ecoli cgmlst schema
ecoli_unpack_cgmlst_schema: $(ECOLI_CGMLST_DIR)/alleles/b4383.fasta

$(ECOLI_CGMLST_DIR)/alleles/b4383.fasta: $(ECOLI_CGMLST_DIR)/alleles/ecoli_cgmlst_alleles_5064703.zip
	$(call log_message,"Unpacking E. coli cgMLST schema ...")
	cd $(ECOLI_CGMLST_DIR)/alleles \
	&& unzip -Duq $$(basename $<) |& tee -a $(INSTALL_LOG)


# Prepping Ecoli cgmlst cgmlst.org schema
ecoli_prep_ecoli_cgmlst_schema: $(ECOLI_CGMLST_DIR)/alleles_rereffed/Escherichia_coli.trn

$(ECOLI_CGMLST_DIR)/alleles_rereffed/Escherichia_coli.trn: $(ECOLI_CGMLST_DIR)/alleles/b4383.fasta
	$(call log_message,"Prepping E. coli cgMLST schema ... WARNING: This takes a looong time. Put on some coffee")
	cd $(ECOLI_CGMLST_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/chewbbaca.sif \
		chewie PrepExternalSchema \
		-i $(ECOLI_CGMLST_DIR)/alleles \
		-o $(ECOLI_CGMLST_DIR)/alleles_rereffed \
		--cpu 1 \
		--ptf $(PRODIGAL_TRAINING_DIR)/Escherichia_coli.trn |& tee -a $(INSTALL_LOG)


# -----------------------------
# K. pneumoniae
# -----------------------------

.PHONY: kpneumoniae_all \
	kpneumoniae_download_reference \
	kpneumoniae_index_reference \
	kpneumoniae_download_cgmlst_schema \
	kpneumoniae_prep_cgmlst_schema

kpneumoniae_all: kpneumoniae_download_reference \
	kpneumoniae_index_reference \
	kpneumoniae_download_cgmlst_schema \
	kpneumoniae_prep_cgmlst_schema


KPNEU_GENOMES_DIR := $(ASSETS_DIR)/genomes/klebsiella_pneumoniae
KPNEU_CGMLST_DIR := $(ASSETS_DIR)/cgmlst/klebsiella_pneumoniae
KPNEU_REFSEQ_ACC := NC_016845.1


kpneumoniae_download_reference: $(KPNEU_GENOMES_DIR)/$(KPNEU_REFSEQ_ACC).fasta

$(KPNEU_GENOMES_DIR)/$(KPNEU_REFSEQ_ACC).fasta:
	$(call log_message,"Downloading K pneumoniae genome ...")
	cd $(SCRIPT_DIR) \
	&& mkdir -p $(KPNEU_GENOMES_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/pythonScripts.sif \
		python3 bin/download_ncbi.py \
		-i $(KPNEU_REFSEQ_ACC) \
		-o $(KPNEU_GENOMES_DIR) |& tee -a $(INSTALL_LOG)


kpneumoniae_index_reference: $(KPNEU_GENOMES_DIR)/$(KPNEU_REFSEQ_ACC).fasta.bwt

$(KPNEU_GENOMES_DIR)/$(KPNEU_REFSEQ_ACC).fasta.bwt: $(KPNEU_GENOMES_DIR)/$(KPNEU_REFSEQ_ACC).fasta
	$(call log_message,"Indexing K pneumoniae genome ...")
	cd $(KPNEU_GENOMES_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/bwakit.sif \
		bwa index $$(basename $<) |& tee -a $(INSTALL_LOG)


# Download Kpneumoniae cgmlst cgmlst.org schema
kpneumoniae_download_cgmlst_schema: $(KPNEU_CGMLST_DIR)/alleles/cgmlst_schema_2187931.zip

$(KPNEU_CGMLST_DIR)/alleles/cgmlst_schema_2187931.zip:
	$(call log_message,"Downloading K. pneumoniae cgMLST schema ...")
	mkdir -p $(KPNEU_CGMLST_DIR)/alleles \
	&& cd $(KPNEU_CGMLST_DIR)/alleles \
	&& wget -O $$(basename $@) https://www.cgmlst.org/ncs/schema/2187931/alleles/ --no-check-certificate |& tee -a $(INSTALL_LOG)


kpneumoniae_unpack_cgmlst_schema: $(KPNEU_CGMLST_DIR)/alleles/KP1_RS24625.fasta

$(KPNEU_CGMLST_DIR)/alleles/KP1_RS24625.fasta: $(KPNEU_CGMLST_DIR)/alleles/cgmlst_schema_2187931.zip
	$(call log_message,"Unpacking K. pneumoniae cgMLST schema ...")
	cd $(KPNEU_CGMLST_DIR)/alleles \
	&& unzip -Duq $$(basename $<) |& tee -a $(INSTALL_LOG)


# Prep Kpneumoniae cgmlst cgmlst.org schema
kpneumoniae_prep_cgmlst_schema: $(KPNEU_CGMLST_DIR)/alleles_rereffed/Klebsiella_pneumoniae.trn

$(KPNEU_CGMLST_DIR)/alleles_rereffed/Klebsiella_pneumoniae.trn:
	$(call log_message,"Prepping K. pneumoniae cgMLST schema ... Warning: This takes a looong time. Put on some coffee!")
	mkdir -p $(KPNEU_CGMLST_DIR)/alleles_rereffed \
	&& cd $(KPNEU_CGMLST_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/chewbbaca.sif \
		chewie PrepExternalSchema \
		-i $(KPNEU_CGMLST_DIR)/alleles \
		-o $(KPNEU_CGMLST_DIR)/alleles_rereffed \
		--cpu 1 \
		--ptf $(PRODIGAL_TRAINING_DIR)/Klebsiella_pneumoniae.trn |& tee -a $(INSTALL_LOG) #TODO: create Kpneumoniae training file

# -----------------------------
# M. tuberculosis
# -----------------------------
MTUBE_GENOMES_DIR := $(ASSETS_DIR)/genomes/mycobacterium_tuberculosis
MTUBE_REFSEQ_ACC := NC_000962.3

.PHONY: mtuberculosis_all \
	mtuberculosis_download_reference\
	mtuberculosis_index_reference

mtuberculosis_all: mtuberculosis_download_reference mtuberculosis_index_reference

mtuberculosis_download_reference: $(MTUBE_GENOMES_DIR)/$(MTUBE_REFSEQ_ACC).fasta

$(MTUBE_GENOMES_DIR)/$(MTUBE_REFSEQ_ACC).fasta:
	$(call log_message,"Downloading M. tuberculosis genome ...")
	mkdir -p $(MTUBE_GENOMES_DIR) \
	&& cd $(SCRIPT_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/pythonScripts.sif \
		python3 bin/download_ncbi.py \
		-i $(MTUBE_REFSEQ_ACC) \
		-o $(MTUBE_GENOMES_DIR) |& tee -a $(INSTALL_LOG)

mtuberculosis_index_reference: $(MTUBE_GENOMES_DIR)/$(MTUBE_REFSEQ_ACC).fasta.bwt

$(MTUBE_GENOMES_DIR)/$(MTUBE_REFSEQ_ACC).fasta.bwt: $(MTUBE_GENOMES_DIR)/$(MTUBE_REFSEQ_ACC).fasta
	$(call log_message,"Indexing M. tuberculosis genome ...")
	cd $(MTUBE_GENOMES_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/bwakit.sif \
		bwa index $$(basename $<) |& tee -a $(INSTALL_LOG)

# ==============================================================================
# Perform checks
# ==============================================================================
# TODO: Integrate with main installation rules?

# -----------------------------
# Check chewBBACA
# -----------------------------
check_chewbbaca:
	@cd $(SCRIPT_DIR) \
	&& saureus=$(SAUR_CGMLST_DIR)/alleles_rereffed \
	&& ecoli=$(ECOLI_CGMLST_DIR)/alleles_rereffed \
	&& kpneumoniae=$(KPNEU_CGMLST_DIR)/alleles_rereffed \
	&& if [[ -d "$$saureus" && -d "$$ecoli" && -d "$$kpneumoniae" ]]; then \
		echo "[✓] PASSED check for chewBBACA: Directories exist:" |& tee -a $(INSTALL_LOG) \
		&& echo "- "$$saureus |& tee -a $(INSTALL_LOG) \
		&& echo "- "$$ecoli |& tee -a $(INSTALL_LOG) \
		&& echo "- $$kpneumoniae" |& tee -a $(INSTALL_LOG); \
	else \
		echo "[!] FAILED check for chewBBACA: Some directories do not exist:"; \
		if [[ ! -d $$saureus ]]; then \
			echo "    Missing directory: $$saureus" |& tee -a $(INSTALL_LOG);  \
		fi; \
		if [[ ! -d $$ecoli ]]; then \
			echo "    Missing directory: $$ecoli" |& tee -a $(INSTALL_LOG);  \
		fi; \
		if [[ ! -d $$kpneumoniae ]]; then \
			echo "    Missing directory: $$kpneumoniae" |& tee -a $(INSTALL_LOG);  \
		fi; \
		echo "    Please report this in an issue on the JASEN repo: https://github.com/genomic-medicine-sweden/JASEN/issues" |& tee -a $(INSTALL_LOG); \
	fi

# -----------------------------
# Check BWA
# -----------------------------
check_bwa:
	@cd $(SCRIPT_DIR) \
	&& ref=$(SAUR_GENOMES_DIR)/$(SAUR_REFSEQ_ACC).fasta; \
		refamb=$${ref}.amb; \
		refann=$${ref}.ann; \
		refbwt=$${ref}.bwt; \
		refpac=$${ref}.pac; \
		refsa=$${ref}.sa \
	&& if [[ -f $$ref \
		&& -f $${refamb} \
		&& -f $${refann} \
		&& -f $${refbwt} \
		&& -f $${refpac} \
		&& -f $${refsa} \
	]]; then \
		echo "[✓] PASSED check for bwa: Indexes exist in $(SAUR_GENOMES_DIR)" |& tee -a $(INSTALL_LOG); \
	else \
		echo "[!] FAILED check for bwa: Indexes do not exist in $(SAUR_GENOMES_DIR)" |& tee -a $(INSTALL_LOG); \
	fi

# -----------------------------
# Check BlastDB
# -----------------------------
MLST_BLAST_DIR := $(ASSETS_DIR)/mlst_db/blast
check_blastdb:
	@cd $(SCRIPT_DIR) \
	&& mlst=$(MLST_BLAST_DIR)/mlst.fa; \
	 mlstndb=$${mlst}.ndb; \
	 mlstnhd=$${mlst}.nhd; \
	 mlstnhi=$${mlst}.nhi; \
	 mlstnhr=$${mlst}.nhr; \
	 mlstnin=$${mlst}.nin; \
	 mlstnog=$${mlst}.nog; \
	 mlstnos=$${mlst}.nos; \
	 mlstnot=$${mlst}.not; \
	 mlstnsq=$${mlst}.nsq; \
	 mlstntf=$${mlst}.ntf; \
	 mlstnto=$${mlst}.nto \
	&& if [[ -f $${mlst} \
		&& -f $${mlstndb} \
		&& -f $${mlstnhd} \
		&& -f $${mlstnhi} \
		&& -f $${mlstnhr} \
		&& -f $${mlstnin} \
		&& -f $${mlstnog} \
		&& -f $${mlstnos} \
		&& -f $${mlstnot} \
		&& -f $${mlstnsq} \
		&& -f $${mlstntf} \
		&& -f $${mlstnto} \
	]]; then \
		echo "[✓] PASSED check for blast: Indexes exist in $(MLST_BLAST_DIR)" |& tee -a $(INSTALL_LOG); \
	else \
		echo "[!] FAILED check for blast: Indexes do not exist in $(MLST_BLAST_DIR)" |& tee -a $(INSTALL_LOG); \
	fi
