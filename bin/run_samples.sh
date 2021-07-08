#!/bin/bash

# http://linuxcommand.org/lc3_man_pages/seth.html
#set -x
# http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
#IFS=$'\n\t'

CWD=$(pwd)
#echo $CWD
SEQ_DATA="$CWD"/"assets/sequencing_data"
#echo $SEQ_DATA
input=$( find assets/sequencing_data/ -name "*.fastq.gz" -and -type f -print0 | xargs -0 echo )
read -r -a input_files_array <<< $input
# Loop through all the isoform csv files
for i in "${!input_files_array[@]}"; do
  FULL_FILE_PATH="$CWD"/$(echo ${input_files_array[$i]})
  RELATIVE_FILE_PATH=$(echo ${input_files_array[$i]})
  #FILE_PATH=$(echo ${input_files_array[$i]} | cut -d "/" -f3-)
  INPUT_DIR=$(echo ${input_files_array[$i]} | awk -F/ '{print $3}')
  FILE_NAME=$(echo ${input_files_array[$i]} | awk -F/ '{print $4}')
  SAMPLE_STEM=$(basename -s '.fastq.gz' $(echo $FILE_NAME))
  SAMPLE_ID=$(echo $SAMPLE_STEM | awk -F_ '{print $1}')
  SPECIES=$(echo "$INPUT_DIR" | awk -F_ '{printf $1"_"$2}')
  FULL_PATH_INPUT_DIR="$SEQ_DATA"/"$INPUT_DIR"

  nextflow run main.nf -profile local,singularity,"$SPECIES" -resume --sample_ID "$INPUT_DIR"
  # Clear reference files before running next workflow
  rm -rf assets/references
  #break
done