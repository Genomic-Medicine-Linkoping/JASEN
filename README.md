# gms-JASEN: Region Östergötland implementation

<!-- TOC -->

- [gms-JASEN: Region Östergötland implementation](#gms-jasen-region-Östergötland-implementation)
    - [Installation](#installation)
        - [Clone and switch to `ro-implementation` branch](#clone-and-switch-to-ro-implementation-branch)
        - [Install Singularity](#install-singularity)
        - [Create Singularity container](#create-singularity-container)
        - [Download reference genomes and create prodigal training files](#download-reference-genomes-and-create-prodigal-training-files)
        - [Move Fastq files and adapter sequences to `assets`](#move-fastq-files-and-adapter-sequences-to-assets)
    - [Usage](#usage)
        - [Optional: Change amount of resources processes are allowed to use](#optional-change-amount-of-resources-processes-are-allowed-to-use)
        - [Run the pipeline](#run-the-pipeline)
        - [Finding results](#finding-results)
- [JASEN](#jasen)
    - [Setup](#setup)
    - [Singularity implementation](#singularity-implementation)
        - [Image creation](#image-creation)
        - [Image execution](#image-execution)
    - [Conda implementation](#conda-implementation)

<!-- /TOC -->

## Installation

### Clone and switch to `ro-implementation` branch

```bash
INSTALL_DIR="dir/where/you/want/your/gms-JASEN/installation"
cd "$INSTALL_DIR"
git clone --recurse-submodules https://github.com/Clinical-Genomics-Linkoping/gms-JASEN.git
cd gms-JASEN
rm -r assets/prodigal_training_files/
git checkout ro-implementation
```

Note: `ro` in the branch name comes from words *Region Östergötland*.

### Install Singularity

Follow the installation instructions: [here](https://sylabs.io/guides/3.7/user-guide/quick_start.html#quick-installation-steps 'Quick installation steps').

If you have RHEL derivative system follow [these instructions](https://sylabs.io/guides/3.0/user-guide/installation.html#install-dependencies 'Installing dependencies with yum/rpm') for installing dependencies.

These instructions were run with globally installed Singularity version 3.7.4.

In case `sudo: singularity: command not found` error is encountered, follow [these instructions](https://sylabs.io/guides/2.5/user-guide/troubleshooting.html#error-running-singularity-with-sudo 'Error running singularity with sudo').

### Create Singularity container

```bash
cd container
bash build_container.sh
```

Note: Building of the container requires sudo privileges.

### Download reference genomes and create prodigal training files

```bash
# Change this line to your specific case
JASEN_INSTALL_DIR="/home/Hanna/Documents/CG-Linkoping/gms-JASEN/"
cd "$JASEN_INSTALL_DIR"
make
```

### Move Fastq files and adapter sequences to `assets`

```bash
cp -r /home/Hanna/Documents/CG-Linkoping/gms-JASEN/assets/sequencing_data/Escherichia_coli_p1 /home/Hanna/Documents/gms-JASEN/assets/sequencing_data/
cp /home/Hanna/Documents/CG-Linkoping/gms-JASEN/assets/adapters/qiaseq_adapters.fa  /home/Hanna/Documents/gms-JASEN/assets/adapters/
```

## Usage

### Optional: Change amount of resources processes are allowed to use

The modifications can be done on lines between 73 and 89 in `nextflow.config`-file.

### Run the pipeline

```bash
make run
```

The pipeline utilises [Resfinder](https://pubmed.ncbi.nlm.nih.gov/22782487/ 'DOI: 10.1093/jac/dks261') to identify acquired antimicrobial resistance genes.


### Finding results

The results can be found in json format in `work/[input-fastq-dir-name]` directory. 

---

**Original instructions for setting up the pipeline:**

---

# JASEN
_Json producing Assembly driven microbial Sequence analysis pipeline to support Epitypification and Normalize classification decisions_

## Setup
* `git clone --recurse-submodules --single-branch --branch master  https://github.com/genomic-medicine-sweden/JASEN.git`
* Edit `JASEN/nextflow.config`
* _`Optionally run: bash JASEN/container/safety_exports.sh USER PREFIX`_


## Singularity implementation
### Image creation
* Install Singularity (through conda or whatever)
* `cd JASEN/container && bash build_container.sh`

### Image execution
* `singularity exec -B JASEN_INSTALL_DIR:/external -B WORKDIR:/out IMAGE nextflow -C /external/nextflow.config run /JASEN/main.nf -profile local,singularity`


## Conda implementation
* Install Conda ( https://www.anaconda.com/distribution )
* Install nextFlow ( `curl -s https://get.nextflow.io | bash` )
* `bash JASEN/setup.sh`
* `nextflow run JASEN/main.nf -profile -local,conda`
