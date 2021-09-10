# JASEN: Region Östergötland implementation

<!-- TOC -->

- [JASEN: Region Östergötland implementation](#jasen-region-östergötland-implementation)
  - [Installation](#installation)
    - [Clone and switch to `ro-implementation` branch](#clone-and-switch-to-ro-implementation-branch)
    - [Install Singularity](#install-singularity)
    - [Download reference genomes and create prodigal training files](#download-reference-genomes-and-create-prodigal-training-files)
    - [Move Fastq.gz files `assets/sequencing_data`](#move-fastqgz-files-assetssequencing_data)
    - [Move adapter sequences to `assets/adapters`](#move-adapter-sequences-to-assetsadapters)
    - [Create a conda environment named `nf`](#create-a-conda-environment-named-nf)
  - [Usage](#usage)
    - [Optional: Change amount of resources Nextflow is allowed to use](#optional-change-amount-of-resources-nextflow-is-allowed-to-use)
    - [Run the pipeline](#run-the-pipeline)
    - [Finding results](#finding-results)
  - [Notes](#notes)

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

### Download reference genomes and create prodigal training files

If you wish to run the pipeline with *Escherichia coli* test data (the fastq input files in: `assets/test_data/sequencing_data/ecoli_1k`), use these commands:

```bash
CONT="library://ljmesi/jasen/main.sif"
make clear_files
cp -r assets/test_data/ref_genomes assets/
cp -r assets/test_data/prodigal_training_files assets/
rm -f assets/prodigal_training_files/Escherichia_coli.trn
cat assets/ref_genomes/Escherichia_coli.fna | singularity exec "$CONT" prodigal -p single -t assets/prodigal_training_files/Escherichia_coli.trn
```

otherwise run this command: 

```bash
make preprocess
```

`make preprocess` removes all downloaded genome files, prodigal training files and checksum file for downloaded genomes (if they previously existed) and then downloads and creates them again. It downloads also the latest spa-typing data from https://spa.ridom.de/

### Move Fastq.gz files `assets/sequencing_data`

The paired end `fastq.gz` files must contain `_R1_` and `_R2_` in the file names in order that Nextflow can recognise the forward and reverse reads.

Here is an example with the *Escherichia coli* test data:

```bash
PROJ_ROOT="$PWD"
INFILES="assets/test_data/sequencing_data/ecoli_1k"
INDIR="assets/sequencing_data/Escherichia_coli_p1/"
cp -r "$PROJ_ROOT"/"$INFILES" "$PROJ_ROOT"/"$INDIR"
```

Here below is an example input file directory for running `make run_samples` command, could be the following:

```
├── sequencing_data
│   ├── Acinetobacter_baumannii_MR130386-1
│   │   ├── MR130386-1_S19_L001_R1_001.fastq.gz
│   │   └── MR130386-1_S19_L001_R2_001.fastq.gz
│   ├── Acinetobacter_baumannii_MR130386-2
│   │   ├── MR130386-2_S20_L001_R1_001.fastq.gz
│   │   └── MR130386-2_S20_L001_R2_001.fastq.gz
│   ├── Acinetobacter_baumannii_NL160346
│   │   ├── NL160346_S17_L001_R1_001.fastq.gz
│   │   └── NL160346_S17_L001_R2_001.fastq.gz
│   ├── Acinetobacter_baumannii_PY1604385
│   │   ├── PY1604385_S18_L001_R1_001.fastq.gz
│   │   └── PY1604385_S18_L001_R2_001.fastq.gz
```

Note 1: All the subdirectories, e.g. `Acinetobacter_baumannii_PY1604385` should be the `sample_ID`s. 

### Move adapter sequences to `assets/adapters`

With the *Escherichia coli* test data we should use *Nextera PE* adapters.

```bash
PROJ_ROOT="$PWD"
TEST_ADAPTERS="assets/test_data/adapters/NexteraPE-PE.fa"
ADAPTERS_DIR="assets/adapters/"
cp "$PROJ_ROOT"/"$TEST_ADAPTERS" "$PROJ_ROOT"/"$ADAPTERS_DIR"
```

### Create a conda environment named `nf`

```bash
conda env create -f nf-env.yml
# or if you have mamba installed:
mamba env create -f nf-env.yml
```

This conda environment contains only the latest version of Nextflow. When the pipeline is run with the `Makefile`, the command preceding the running of the pipeline, is activating this environment.  

### Create Kaiju and Kraken2 databases

The databases can be downloaded with and built with commands:
```bash
KAI_DB=[Absolute-path-to-db]
KRA_DB=[Absolute-path-to-dir-where-the-db-dir-named-KRA_DB_NAME-will-be-created]
KRA_DB_NAME=K2DB
make create_kaijudb KAIJU_DB="$KAI_DB"
make create_kraken2db \
KRAKEN_DB_DIR="$KRA_DB" \
KRAKEN_DB_NAME="$KRA_DB_NAME"
```

## Usage

### Optional: Change amount of resources Nextflow is allowed to use

Use the following command for finding which lines to modify in order to adjust (in the `nextflow.config`-file) the resource usage according to your local system:

```bash
grep -nA 16 -P "^process\s\{" nextflow.config
```

### Run the pipeline

The pipeline can be run with the following command (after having performed the steps in [Installation](#installation)).

Adjust the values for `SPECIES` and `SAMPLE_ID` according to your case. The values should be:

- `SPECIES` = The name of the species the fastq samples are from
- `SAMPLE_ID` = The input directory name inside `assets/sequencing_data`

Below is one example of running the pipeline with test data

```bash
make \
SPECIES=Escherichia_coli \
SAMPLE_ID=Escherichia_coli_p1 \
ARGS="--genome_name Escherichia_coli --adapter_fname NexteraPE-PE.fa --prodigal_file Escherichia_coli"
```

Note 1: The underscores between words in make commandline arguments `SPECIES` and `SAMPLE_ID` are mandatory.

Note 2: The last line (starting with `ARGS=`) is useful when testing the pipeline. It is not necessary otherwise. 

Note 3: `prodigal_file` is the just the basename (i.e. without `.trn`) and should be located in `assets/prodigal_training_files/`, e.g. `assets/prodigal_training_files/Escherichia_coli.trn`.

### Finding results

The results can be found in json format in `results/[input-fastq-dir-name]` directory. 

## Notes

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

