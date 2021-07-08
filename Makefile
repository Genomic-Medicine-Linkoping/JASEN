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

download_bacterial_genomes:
	@echo ""
	@echo "Laddar ner $(MT) refseq genomet, se mera info om genomet här: https://www.ncbi.nlm.nih.gov/assembly/GCF_000195955.2/"
	@echo ""
	wget -vc -O $(RG)/Mycobacterium_tuberculosis_GCF_000195955.2_ASM19595v2.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Mycobacterium_tuberculosis/reference/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz
	@echo ""
	@echo "Skapar md5sum av $(MT) refseq genomet:"
	md5sum $(RG)/Mycobacterium_tuberculosis_GCF_000195955.2_ASM19595v2.fna.gz > $(RG)/md5sums.txt
	@echo ""
	@echo ""
	@echo "Laddar ner $(KP) refseq genomet, se mera info om genomet här: https://www.ncbi.nlm.nih.gov/assembly/GCF_000240185.1/"
	@echo ""
	wget -vc -O $(RG)/Klebsiella_pneumoniae_GCF_000240185.1_ASM24018v2.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Klebsiella_pneumoniae/reference/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna.gz
	@echo ""
	@echo "Skapar md5sum av $(KP) refseq genomet:"
	md5sum $(RG)/Klebsiella_pneumoniae_GCF_000240185.1_ASM24018v2.fna.gz >> $(RG)/md5sums.txt
	@echo ""
	@echo ""
	@echo "Laddar ner $(SA) refseq genomet, se mera info om genomet här: https://www.ncbi.nlm.nih.gov/assembly/GCF_000013425.1/"
	wget -vc -O $(RG)/Staphylococcus_aureus_GCF_000013425.1_ASM1342v1.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Staphylococcus_aureus/reference/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz
	@echo ""
	@echo "Skapar md5sum av $(SA) refseq genomet"
	@echo ""
	md5sum $(RG)/Staphylococcus_aureus_GCF_000013425.1_ASM1342v1.fna.gz >> $(RG)/md5sums.txt
	@echo ""
	@echo ""
	@echo "Laddar ner $(EC) refseq genomet, se mera info om genomet här: https://www.ncbi.nlm.nih.gov/assembly/GCF_000008865.2/"
	wget -vc -O $(RG)/Escherichia_coli_GCF_000008865.2_ASM886v2.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/reference/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz
	@echo ""
	@echo "Skapar md5sum av $(EC) refseq genomet"
	@echo ""
	md5sum $(RG)/Escherichia_coli_GCF_000008865.2_ASM886v2.fna.gz >> $(RG)/md5sums.txt
	@echo ""
	@echo ""
	@echo "Laddar ner $(AC) refseq genomet, se mera info om genomet här: https://www.ncbi.nlm.nih.gov/genome/403"
	wget -vc -O $(RG)/Acinetobacter_baumannii_GCF_002116925.1_ASM211692v1.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/116/925/GCF_002116925.1_ASM211692v1/GCF_002116925.1_ASM211692v1_genomic.fna.gz
	@echo ""
	@echo "Skapar md5sum av $(AC) refseq genomet"
	@echo ""
	md5sum $(RG)/Acinetobacter_baumannii_GCF_002116925.1_ASM211692v1.fna.gz >> $(RG)/md5sums.txt
	@echo ""
	@echo ""
	@echo "Laddar ner $(ES) refseq genomet, se mera info om genomet här: https://www.ncbi.nlm.nih.gov/genome/808?genome_assembly_id=389507"
	wget -vc -O $(RG)/Enterococcus_faecalis_GCF_003319815.1_ASM331981v1.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/319/815/GCF_003319815.1_ASM331981v1/GCF_003319815.1_ASM331981v1_genomic.fna.gz
	@echo ""
	@echo "Skapar md5sum av $(ES) refseq genomet"
	@echo ""
	md5sum $(RG)/Enterococcus_faecalis_GCF_003319815.1_ASM331981v1.fna.gz >> $(RG)/md5sums.txt
	@echo ""
	@echo ""
	@echo "Laddar ner $(EF) refseq genomet, se mera info om genomet här: https://www.ncbi.nlm.nih.gov/genome/871?genome_assembly_id=429798"
	wget -vc -O $(RG)/Enterococcus_faecium_GCF_000336405.1_ASM33640v1.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/336/405/GCF_000336405.1_ASM33640v1/GCF_000336405.1_ASM33640v1_genomic.fna.gz
	@echo ""
	@echo "Skapar md5sum av $(EF) refseq genomet"
	@echo ""
	md5sum $(RG)/Enterococcus_faecium_GCF_000336405.1_ASM33640v1.fna.gz >> $(RG)/md5sums.txt
	@echo ""
	@echo ""
	@echo "Laddar ner $(PA) refseq genomet, se mera info om genomet här: https://www.ncbi.nlm.nih.gov/genome/187?genome_assembly_id=299953"
	wget -vc -O $(RG)/Pseudomonas_aeruginosa_GCF_000006765.1_ASM676v1.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz
	@echo ""
	@echo "Skapar md5sum av $(PA) refseq genomet"
	@echo ""
	md5sum $(RG)/Pseudomonas_aeruginosa_GCF_000006765.1_ASM676v1.fna.gz >> $(RG)/md5sums.txt
	@echo ""
	@echo ""
	@echo "Laddar ner $(CD) refseq genomet, se mera info om genomet här: https://www.ncbi.nlm.nih.gov/genome/535?genome_assembly_id=308203"
	wget -vc -O $(RG)/Clostridioides_difficile_GCF_002007885.1_ASM200788v1.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/007/885/GCF_002007885.1_ASM200788v1/GCF_002007885.1_ASM200788v1_genomic.fna.gz
	@echo ""
	@echo "Skapar md5sum av $(CD) refseq genomet"
	@echo ""
	md5sum $(RG)/Clostridioides_difficile_GCF_002007885.1_ASM200788v1.fna.gz >> $(RG)/md5sums.txt
	@echo ""
	@echo ""
	@echo "Laddar ner $(MA) refseq genomet, se mera info om genomet här: https://www.ncbi.nlm.nih.gov/genome/?term=Mycobacterium+africanum"
	wget -vc -O $(RG)/Mycobacterium_africanum_GCF_000253355.1_ASM25335v1.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/253/355/GCF_000253355.1_ASM25335v1/GCF_000253355.1_ASM25335v1_genomic.fna.gz
	@echo ""
	@echo "Skapar md5sum av $(MA) refseq genomet"
	@echo ""
	md5sum $(RG)/Mycobacterium_africanum_GCF_000253355.1_ASM25335v1.fna.gz >> $(RG)/md5sums.txt
	@echo ""
	@echo ""
	@echo "Laddar ner $(MB) refseq genomet, se mera info om genomet här: https://www.ncbi.nlm.nih.gov/genome/?term=Mycobacterium+bovis"
	wget -vc -O $(RG)/Mycobacterium_bovis_GCF_005156105.1_ASM515610v1.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/005/156/105/GCF_005156105.1_ASM515610v1/GCF_005156105.1_ASM515610v1_genomic.fna.gz
	@echo ""
	@echo "Skapar md5sum av $(MB) refseq genomet"
	@echo ""
	md5sum $(RG)/Mycobacterium_bovis_GCF_005156105.1_ASM515610v1.fna.gz >> $(RG)/md5sums.txt
	@echo ""
	@echo ""
	@echo "Laddar ner $(SE) refseq genomet, se mera info om genomet här: https://www.ncbi.nlm.nih.gov/genome/152?genome_assembly_id=305903"
	wget -vc -O $(RG)/Salmonella_enterica_GCF_000006945.2_ASM694v2.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz
	@echo ""
	@echo "Skapar md5sum av $(SE) refseq genomet"
	@echo ""
	md5sum $(RG)/Salmonella_enterica_GCF_000006945.2_ASM694v2.fna.gz >> $(RG)/md5sums.txt
	@echo ""

create_prodigal_trn_files:
	@echo ""
	@echo "Skapar mapp för singularity körningen:"
	@mkdir -p work
	@echo ""
	@echo "Skapar prodigal träningsfiler från refseq genomen:"
	@echo ""
	@echo "För $(MT):"
	zcat $(RG)/Mycobacterium_tuberculosis_GCF_000195955.2_ASM19595v2.fna.gz | $(SG) $(PT)/Mycobacterium_tuberculosis_GCF_000195955.2_ASM19595v2.trn
	@echo ""
	@echo ""
	@echo "För $(KP):"
	zcat $(RG)/Klebsiella_pneumoniae_GCF_000240185.1_ASM24018v2.fna.gz | $(SG) $(PT)/Klebsiella_pneumoniae_GCF_000240185.1_ASM24018v2.trn
	@echo ""
	@echo ""
	@echo "För $(SA):"
	zcat $(RG)/Staphylococcus_aureus_GCF_000013425.1_ASM1342v1.fna.gz | $(SG) $(PT)/Staphylococcus_aureus_GCF_000013425.1_ASM1342v1.trn
	@echo ""
	@echo ""
	@echo "För $(EC):"
	zcat $(RG)/Escherichia_coli_GCF_000008865.2_ASM886v2.fna.gz | $(SG) $(PT)/Escherichia_coli_GCF_000008865.2_ASM886v2.trn
	@echo ""
	@echo ""
	@echo "För $(AC):"
	zcat $(RG)/Acinetobacter_baumannii_GCF_002116925.1_ASM211692v1.fna.gz | $(SG) $(PT)/Acinetobacter_baumannii_GCF_002116925.1_ASM211692v1.trn
	@echo ""
	@echo ""
	@echo "För $(ES):"
	zcat $(RG)/Enterococcus_faecalis_GCF_003319815.1_ASM331981v1.fna.gz | $(SG) $(PT)/Enterococcus_faecalis_GCF_003319815.1_ASM331981v1.trn
	@echo ""
	@echo ""
	@echo "För $(EF):"
	zcat $(RG)/Enterococcus_faecium_GCF_000336405.1_ASM33640v1.fna.gz | $(SG) $(PT)/Enterococcus_faecium_GCF_000336405.1_ASM33640v1.trn
	@echo ""
	@echo ""
	@echo "För $(PA):"
	zcat $(RG)/Pseudomonas_aeruginosa_GCF_000006765.1_ASM676v1.fna.gz | $(SG) $(PT)/Pseudomonas_aeruginosa_GCF_000006765.1_ASM676v1.trn
	@echo ""
	@echo ""
	@echo "För $(CD):"
	zcat $(RG)/Clostridioides_difficile_GCF_002007885.1_ASM200788v1.fna.gz | $(SG) $(PT)/Clostridioides_difficile_GCF_002007885.1_ASM200788v1.trn
	@echo ""
	@echo ""
	@echo "För $(MA):"
	zcat $(RG)/Mycobacterium_africanum_GCF_000253355.1_ASM25335v1.fna.gz | $(SG) $(PT)/Mycobacterium_africanum_GCF_000253355.1_ASM25335v1.trn
	@echo ""
	@echo ""
	@echo "För $(MB):"
	zcat $(RG)/Mycobacterium_bovis_GCF_005156105.1_ASM515610v1.fna.gz | $(SG) $(PT)/Mycobacterium_bovis_GCF_005156105.1_ASM515610v1.trn
	@echo ""
	@echo ""
	@echo "För $(SE):"
	zcat $(RG)/Salmonella_enterica_GCF_000006945.2_ASM694v2.fna.gz | $(SG) $(PT)/Salmonella_enterica_GCF_000006945.2_ASM694v2.trn
	@echo ""

uncompress_genomes:
	zcat $(RG)/Mycobacterium_tuberculosis_GCF_000195955.2_ASM19595v2.fna.gz > $(RG)/Mycobacterium_tuberculosis_GCF_000195955.2_ASM19595v2.fna
	zcat $(RG)/Klebsiella_pneumoniae_GCF_000240185.1_ASM24018v2.fna.gz > $(RG)/Klebsiella_pneumoniae_GCF_000240185.1_ASM24018v2.fna
	zcat $(RG)/Staphylococcus_aureus_GCF_000013425.1_ASM1342v1.fna.gz > $(RG)/Staphylococcus_aureus_GCF_000013425.1_ASM1342v1.fna
	zcat $(RG)/Escherichia_coli_GCF_000008865.2_ASM886v2.fna.gz > $(RG)/Escherichia_coli_GCF_000008865.2_ASM886v2.fna
	zcat $(RG)/Acinetobacter_baumannii_GCF_002116925.1_ASM211692v1.fna.gz > $(RG)/Acinetobacter_baumannii_GCF_002116925.1_ASM211692v1.fna
	zcat $(RG)/Enterococcus_faecalis_GCF_003319815.1_ASM331981v1.fna.gz > $(RG)/Enterococcus_faecalis_GCF_003319815.1_ASM331981v1.fna
	zcat $(RG)/Enterococcus_faecium_GCF_000336405.1_ASM33640v1.fna.gz > $(RG)/Enterococcus_faecium_GCF_000336405.1_ASM33640v1.fna
	zcat $(RG)/Pseudomonas_aeruginosa_GCF_000006765.1_ASM676v1.fna.gz > $(RG)/Pseudomonas_aeruginosa_GCF_000006765.1_ASM676v1.fna
	zcat $(RG)/Clostridioides_difficile_GCF_002007885.1_ASM200788v1.fna.gz > $(RG)/Clostridioides_difficile_GCF_002007885.1_ASM200788v1.fna
	zcat $(RG)/Mycobacterium_africanum_GCF_000253355.1_ASM25335v1.fna.gz > $(RG)/Mycobacterium_africanum_GCF_000253355.1_ASM25335v1.fna
	zcat $(RG)/Mycobacterium_bovis_GCF_005156105.1_ASM515610v1.fna.gz > $(RG)/Mycobacterium_bovis_GCF_005156105.1_ASM515610v1.fna
	zcat $(RG)/Salmonella_enterica_GCF_000006945.2_ASM694v2.fna.gz > $(RG)/Salmonella_enterica_GCF_000006945.2_ASM694v2.fna

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
	# TODO: Could add the rest of the steps here as well
	# /usr/bin/git add assets/var-genes-ro
	# /usr/bin/git push $(UPSTR_NAME) $(CURR_BRANCH)
