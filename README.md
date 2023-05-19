<p align="center">
  <a href="https://github.com/genomic-medicine-sweden/JASEN">
    <img src="artwork/logo.png"/>
  </a>
</p>

_Json producing Assembly driven microbial Sequence analysis pipeline to support Epitypification and Normalize classification decisions_

JASEN produces results for epidemiological and surveillance purposes.
JASEN has been tested using MRSA, but should work well with any bacteria with a stable cgMLST scheme.

## Requirements

* Singularity
* Nextflow (`curl -s https://get.nextflow.io | bash`)

### Recommended

* Conda
* Singularity Remote Login

## Development deployment (self-contained)

### Copy code locally

```
git clone --recurse-submodules --single-branch --branch master  https://github.com/genomic-medicine-sweden/JASEN.git && cd JASEN
```

### Create conda environment needed for `deploy_references.sh`

```
bash -i deploy/deploy_conda.sh
```

### Download references and databases. NOTE: Ensure that after running `deploy_references.sh` that there are no error messages

```
bash -i deploy/deploy_references.sh
```

### Download references and databases using singularity instead of conda (above). NOTE: Ensure that after running `deploy_references_singularity.sh` that there are no error messages

```
bash -i deploy/deploy_references_singularity.sh
```

### Access to OCI regestries (Optional)

```
singularity remote login
```

### Creates singularity images. NOTE: Ensure that you have sudo priviledges before running `build_container.sh`

```
cd container && sudo bash -i build_container.sh && cd ..
```

## Config and test data

### Config (`configs/nextflow.base.config`)

* Edit the `root` parameter in `configs/nextflow.base.config`
* Edit the `krakenDb`, `workDir` and `outdir` parameters in `configs/nextflow.base.config`
* Edit the `runOptions` in `configs/nextflow.base.config` in order to mount directories to your run

### Test data (`assets/test_data/samplelist.csv`)

* Edit the read1 and read2 columns in `assets/test_data/samplelist.csv`

## Setting temp directories

### Open `~/.bashrc` and add the following (Edit `/path/to/tmp`)

```
export SINGULARITY_TMPDIR="/path/to/tmp"
```

## Fetching databases

### Choose between MiniKraken DB (8GB) or Kraken DB (64GB; Recommended***)

### Download MiniKraken database

```
wget -O /path/to/kraken_db/krakenmini.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20230314.tar.gz
tar -xf /path/to/kraken_db/krakenmini.tar.gz
```

### Download Kraken database

```
wget -O /path/to/kraken_db/krakenstd.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230314.tar.gz
tar -xf /path/to/kraken_db/krakenstd.tar.gz
```

## Usage

### Simple self-test

```
nextflow run main.nf -entry bacterial_default -profile staphylococcus_aureus -config configs/nextflow.base.config --csv assets/test_data/samplelist.csv
```

### Usage arguments

| Argument type | Options                                | Required |
| ------------- | -------------------------------------- | -------- |
| -profile      | staphylococcus_aureus/escherichia_coli | True     |
| -entry        | bacterial_default                      | True     |
| -config       | nextflow.base.config                   | True     |
| -resume       | NA                                     | False    |
| --output      | user specified                         | False    |

### Input file format 

```csv
id,species,platform,read1,read2
p1,saureus,illumina,assets/test_data/sequencing_data/saureus_10k/saureus_large_R1_001.fastq.gz,assets/test_data/sequencing_data/saureus_10k/saureus_large_R2_001.fastq.gz
```

### Input species options

* `saureus`
* `ecoli`

## Component Breakdown

### QC

* [Kraken2](https://ccb.jhu.edu/software/kraken2/): Species detection.
* [Bracken](https://ccb.jhu.edu/software/bracken/): Combined with Kraken2 for species detection.
* [bwa mem](https://github.com/lh3/bwa): Maps reads to cgMLST loci (demarcated by bed file) in order to estimate genome coverage. Low levels of Intra-species contamination or erroneous mapping is removed using bwa and filtering away the heterozygous mapped bases.
* [interquartile range](https://en.wikipedia.org/wiki/Interquartile_range): Calculates evenness of coverage.

### Assembly

* [SPAdes](http://cab.spbu.ru/software/spades/): De novo assembly for Ion Torrent.
* [SKESA](https://www.ridom.de/seqsphere/ug/v60/SKESA_Assembler.html): De novo assembly for Illumina.
* [QUAST](http://cab.spbu.ru/software/quast/): Extracts QC data from the assembly.

### Epidemiological typing

* [chewBBACA](https://github.com/B-UMMI/chewBBACA/wiki): Calculates cgMLST of extracted alleles decided by schema. Number of missing loci is calculated and used as a QC parameter.
* [cgmlst.net](https://www.cgmlst.org/ncs/schema/141106/): The cgMLST reference schema.
* [mlst](https://github.com/tseemann/mlst): Caculates traditional 7-locus MLST.

#### Currently, 2 profiles are supported:

* `staphylococcus_aureus` (`saureus` in csv)
* `escherichia_coli` (`ecoli` in csv)

#### Future profiles that will be supported:

* `klebsiella_pneumoniae`
* `mycobacterium_tuberculosis`

### Virulence and resistance markers

* [resfinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/): Detects antimicrobial resistance genes as well as environmental and chemical resistance genes.
* [pointfinder](https://bitbucket.org/genomicepidemiology/pointfinder/src/master/): Combines with resfinder to detect variants.
* [virulencefinder](https://bitbucket.org/genomicepidemiology/virulencefinder/src/master/): Detects virulence genes.
* [amrfinderplus](https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus): Detects antimicrobial resistance genes as well as environmental, chemical resistance and virulence genes.
* [resfinder_db](https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/): Resfinder database.
* [pointfinder_db](https://bitbucket.org/genomicepidemiology/pointfinder_db/src/master/): Pointfinder database.
* [virulencefinder_db](https://bitbucket.org/genomicepidemiology/virulencefinder_db/src/master/): Virulencefinder database.

## Report and visualisation

* [Bonsai](https://github.com/Clinical-Genomics-Lund/cgviz): Visualises JASEN outputs.
* [graptetree](https://github.com/achtman-lab/GrapeTree): Visualise phylogenetic relationship using cgmlst data.

## Tips

It is recommended that you use latest versions of software tools, however if you are running an older version of Singularity and you get an error `FATAL: could not open image JASEN/container/*.sif: image format not recognized!` check the permissions set on image `*.sif`. Make sure you have the permission to execute it.
