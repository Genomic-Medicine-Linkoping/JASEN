#!/bin/bash

# http://linuxcommand.org/lc3_man_pages/seth.html
#set -x
# http://redsymbol.net/articles/unofficial-bash-strict-mode/
#set -e
set -uo pipefail
#IFS=$'\n\t'

CWD=$(pwd)
SEQ_DATA="$CWD"/"assets/sequencing_data"

# Clear log file for the total run
rm -f nf_logs.log

input=$( find assets/sequencing_data/ -name "*R1_001.fastq.gz" -and -type f -print0 | xargs -0 echo )
read -r -a input_files_array <<< $input
# Loop through all the samples
for i in "${!input_files_array[@]}"; do
  # /home/rada/Documents/CGL/JASEN/assets/sequencing_data/Proteus_vulgaris_MRR160070/MRR160070_S14_L001_R1_001.fastq.gz
  FULL_FILE_PATH="$CWD"/$(echo ${input_files_array[$i]})
  # assets/sequencing_data/Proteus_vulgaris_MRR160070/MRR160070_S14_L001_R1_001.fastq.gz
  RELATIVE_FILE_PATH=$(echo ${input_files_array[$i]})
  # Proteus_vulgaris_MRR160070
  INPUT_DIR=$(echo ${input_files_array[$i]} | awk -F/ '{print $3}')
  # MRR160070_S14_L001_R1_001.fastq.gz
  FILE_NAME=$(echo ${input_files_array[$i]} | awk -F/ '{print $4}')
  # MRR160070_S14_L001_R1_001
  SAMPLE_STEM=$(basename -s '.fastq.gz' $(echo $FILE_NAME))
  # MRR160070
  SAMPLE_ID=$(echo $SAMPLE_STEM | awk -F_ '{print $1}')
  # Proteus_vulgaris
  SPECIES=$(echo "$INPUT_DIR" | awk -F_ '{printf $1"_"$2}')
  # /home/rada/Documents/CGL/JASEN/assets/sequencing_data/Proteus_vulgaris_MRR160070
  FULL_PATH_INPUT_DIR="$SEQ_DATA"/"$INPUT_DIR"
  # Clear reference files before running next workflow run
  rm -rf assets/references
  start_time=$(date +"%c")
  echo "Start time: $start_time" | tee -a nf_logs.log
  start=$(date +%s)
  nextflow run main.nf -profile local,singularity,"$SPECIES" --sample_ID "$INPUT_DIR" | tee -a nf_logs.log
  end=$(date +%s)
  runtime=$((end-start))
  echo "The run took: $runtime s." | tee -a nf_logs.log
  end_time=$(date +"%c")
  echo "End time: $end_time" | tee -a nf_logs.log
  #break
done