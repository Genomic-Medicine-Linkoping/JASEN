# JASEN: Region Östergötland implementation

<!-- TOC -->

- [JASEN: Region Östergötland implementation](#jasen-region-östergötland-implementation)
  - [Installation](#installation)
    - [Clone and switch to `ro-implementation` branch](#clone-and-switch-to-ro-implementation-branch)
    - [Install Singularity](#install-singularity)
    - [Create Singularity containers](#create-singularity-containers)
    - [Download reference genomes and create prodigal training files](#download-reference-genomes-and-create-prodigal-training-files)
    - [Move Fastq.gz files `assets/sequencing_data`](#move-fastqgz-files-assetssequencing_data)
    - [Move adapter sequences to `assets/adapters`](#move-adapter-sequences-to-assetsadapters)
    - [Create a conda environment named `nf`](#create-a-conda-environment-named-nf)
    - [Notes](#notes)
  - [Usage](#usage)
    - [Optional: Change amount of resources Nextflow is allowed to use](#optional-change-amount-of-resources-nextflow-is-allowed-to-use)
    - [Run the pipeline](#run-the-pipeline)
    - [Finding results](#finding-results)

<!-- /TOC -->

## Installation

### Clone and switch to `ro-implementation` branch

```bash
INSTALL_DIR="dir/where/you/want/your/gms-JASEN/installation"
cd "$INSTALL_DIR"
git clone --recurse-submodules --single-branch --branch ro-implementation https://github.com/Genomic-Medicine-Linkoping/JASEN.git
cd gms-JASEN
rm -r assets/prodigal_training_files/
git checkout ro-implementation
```

Note: `ro` in the branch name comes from words *Region Östergötland*.

### Install Singularity

Follow the installation instructions: [here](https://sylabs.io/guides/3.8/user-guide/quick_start.html 'Quick installation steps').

If you have RHEL derivative system follow [these instructions](https://sylabs.io/guides/3.0/user-guide/installation.html#install-dependencies 'Installing dependencies with yum/rpm') for installing dependencies.

These instructions were tested with globally installed Singularity version 3.7.4.

In case `sudo: singularity: command not found` error is encountered, follow [these instructions](https://sylabs.io/guides/2.5/user-guide/troubleshooting.html#error-running-singularity-with-sudo 'Error running singularity with sudo').

### Create Singularity containers

```bash
make build_containers
```

Note: Building of the container requires sudo privileges.

This command creates two singularity images which are used in the pipeline. They are:

- `jasen_<date>.sif`: This is used in all processes except in `build_report` (building the final html report)
- `jasen_tidyverse_<date>.sif`: This is used in building the final html report with `build_report` process

### Download reference genomes and create prodigal training files

```bash
make
```

This command 

### Move Fastq.gz files `assets/sequencing_data`

The paired end `fastq.gz` files must contain `_R1_` and `_R2_` in the file names in order that Nextflow can recognise the forward and reverse reads.

```bash
# Fastq-files
# Changes these to your local settings
PROJ_ROOT="/home/rada/Documents/CGL/JASEN"
INFILES="assets/test_data/sequencing_data/ecoli_1k"
INDIR="assets/sequencing_data/"
cp -r "$PROJ_ROOT"/"$INFILES" "$PROJ_ROOT"/"$INDIR"
```

### Move adapter sequences to `assets/adapters`

```bash

TEST_ADAPTERS="assets/test_data/adapters/giaseq_adapters.fa"
ADAPTERS_DIR="assets/adapters/"
cp "$PROJ_ROOT"/"$TEST_ADAPTERS" "$PROJ_ROOT"/"$ADAPTERS_DIR"
```

### Create a conda environment named `nf`

```bash
conda env create -f nf-env.yml
```

The pipeline utilises [Resfinder](https://pubmed.ncbi.nlm.nih.gov/22782487/ 'DOI: 10.1093/jac/dks261') to identify acquired antimicrobial resistance genes.


### Finding results

The results can be found in json format in `work/[input-fastq-dir-name]` directory. 

---

**Original instructions for setting up the pipeline:**

---

# JASEN <!-- omit in toc -->
_Json producing Assembly driven microbial Sequence analysis pipeline to support Epitypification and Normalize classification decisions_

## Setup <!-- omit in toc -->
* `git clone --recurse-submodules --single-branch --branch master  https://github.com/genomic-medicine-sweden/JASEN.git`
* Edit `JASEN/nextflow.config`
* _`Optionally run: bash JASEN/container/safety_exports.sh USER PREFIX`_


## Singularity implementation <!-- omit in toc -->
### Image creation <!-- omit in toc -->
* Install Singularity (through conda or whatever)
* `cd JASEN/container && bash build_container.sh`

### Image execution <!-- omit in toc -->
* `singularity exec -B JASEN_INSTALL_DIR:/external -B WORKDIR:/out IMAGE nextflow -C /external/nextflow.config run /JASEN/main.nf -profile local,singularity`


## Conda implementation <!-- omit in toc -->
* Install Conda ( https://www.anaconda.com/distribution )
* Install nextFlow ( `curl -s https://get.nextflow.io | bash` )
* `bash JASEN/setup.sh`
* `nextflow run JASEN/main.nf -profile -local,conda`
