# JASEN: Region Östergötland implementation

<!-- TOC -->

- [JASEN: Region Östergötland implementation](#jasen-region-östergötland-implementation)
  - [Installation](#installation)
    - [Clone and switch to `ro-implementation` branch](#clone-and-switch-to-ro-implementation-branch)
    - [Install Singularity](#install-singularity)
    - [Build the main Singularity image](#build-the-main-singularity-image)
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
INSTALL_DIR="dir/where/you/want/your/JASEN/installation"
cd "$INSTALL_DIR"
git clone --recurse-submodules --single-branch --branch ro-implementation https://github.com/Genomic-Medicine-Linkoping/JASEN.git
cd JASEN
git checkout ro-implementation
```

Note: `ro` in the branch name comes from words *[Region Östergötland](https://www.regionostergotland.se/)*.

### Install Singularity

Follow the installation instructions: [here](https://sylabs.io/guides/3.8/user-guide/quick_start.html 'Quick installation steps').

If you have RHEL derivative system follow [these instructions](https://sylabs.io/guides/3.0/user-guide/installation.html#install-dependencies 'Installing dependencies with yum/rpm') for installing dependencies.

These instructions were tested with globally installed Singularity version 3.7.4.

In case `sudo: singularity: command not found` error is encountered, follow [these instructions](https://sylabs.io/guides/2.5/user-guide/troubleshooting.html#error-running-singularity-with-sudo 'Error running singularity with sudo').

### Build the main Singularity image

```bash
make build_sif
```

Note: Building of the image requires sudo privileges.

This command creates one singularity image which is used in the pipeline. It is `jasen_<date>.sif`. It is used in all processes except in `build_report` (building the final html report). For creating a final report is used another singularity [image](https://cloud.sylabs.io/library/ljmesi/default/jasen_tidyverse.sif) hosted on sylabs cloud library. `Singularity_tidyverse` image in `container/` directory was used to build this image.

### Download reference genomes and create prodigal training files

```bash
SIF=jasen_<date>.sif
make preprocess CONT_NAME="$SIF"
```

where `jasen_<date>.sif` is the name of the previous step newly built singularity image, e.g. `jasen_2021-07-28.sif`.

This command removes all downloaded genome files, prodigal training files and checksum file (if they previously existed) and then downloads creates them again.

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

### Notes

Note that currently only these species are supported for determining both CGMLST and MLST:
- Escherichia coli
- Staphylococcus aureus
- Klebsiella pneumoniae
- Mycobacterium tuberculosis
- Acinetobacter baumannii
- Enterococcus faecalis
- Enterococcus faecium
- Pseudomonas aeruginosa
- Clostridioides difficile
- Mycobacterium africanum
- Mycobacterium bovis
- Salmonella enterica

and these are supported only for determining MLST:
- Enterobacter cloacae
- Proteus mirabilis
- Proteus vulgaris
- Mycobacterium gordonae
- Mycobacteroides abscessus
- Mycobacterium avium
- Mycobacterium intracellulare
- Mycobacterium malmoense
- Staphylococcus argenteus
- Staphylococcus saprophyticus
- Stenotrophomonas maltophilia
- Streptococcus pyogenes
- Citrobacter braakii
- Corynebacterium striatum
- Enterococcus gallinarum
- Klebsiella aerogenes
- Mycobacterium chimera
- Mycobacterium malmoense
- Mycobacterium kansasii
- Mycobacteroides chelonae
- Mycobacterium celatum
- Mycobacterium marinum
- Mycobacterium szulgai
- Mycobacterium scrofulaceum
- Mycobacterium xenopi


## Usage

### Optional: Change amount of resources Nextflow is allowed to use

Use the following command for finding which lines to modify in order to adjust (in the `nextflow.config`-file) the resource usage according to your local system:

```bash
grep -nA 16 -P "^process\s\{" nextflow.config
```

### Run the pipeline

The pipeline can be run with the following command (after having performed the steps in [Installation](#installation)).

Adjust the values for `SPECIES`, `SAMPLE_ID`, `CONT_NAME` and `CONT_REPORT` according to your case. The values should be:

- `SPECIES` = The name of the species the fastq samples are from
- `SAMPLE_ID` = The input directory name inside `assets/sequencing_data`
- `CONT_NAME` = The name of the singularity image (which does not contain the `tidyverse` word in it) created in [Create Singularity containers](#create-singularity-containers)
- `CONT_REPORT` = The name of the singularity image (which contains the `tidyverse` word in it) created in [Create Singularity containers](#create-singularity-containers)

Below is one example of running the pipeline.  

```bash
make run \
SPECIES=Staphylococcus_saprophyticus \
SAMPLE_ID=Staphylococcus_saprophyticus_Stam-121 \
CONT_NAME=jasen_2021-07-27.sif \
CONT_REPORT=jasen_tidyverse_2021-07-27.sif
```

Note: The underscores between words in make commandline arguments `SPECIES` and `SAMPLE_ID` are mandatory.

### Finding results

The results can be found in json format in `results/[input-fastq-dir-name]` directory. 
