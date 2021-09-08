#!/bin/bash
set -e
set -u
set -o pipefail

# Commands run
KAIJU_DB=/data/CGL/JASEN/kaiju-db
cd "$KAIJU_DB"
kaiju-makedb -t 50 -s nr_euk
mkdir reads out

READS=/home/rada/Documents/TNscope/exp/test-kaiju/reads

mv trim_front_pair.fastqc.gz trim_rev_pair.fastq.gz reads/

kaiju -z 25 -t nodes.dmp -f nr_euk/kaiju_db_nr_euk.fmi -i reads/trim_front_pair.fastq.gz -j reads/trim_rev_pair.fastq.gz -o out/kaiju.out

kaiju2table -t nodes.dmp -n names.dmp -r species -o out/kaiju_summary.tsv out/kaiju.out -v
multiqc out/