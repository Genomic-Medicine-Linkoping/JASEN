#!/bin/bash
set -e
set -u
set -o pipefail

INPUT_FILE_DIR=$1
FIND_PATTERNS=$2
SAMPLE_DATA_IN_FULL_PATH_LOC=$3
OUT_DIR=$4
FILE_TYPE=$5
# echo $INPUT_FILE_DIR
# echo $FIND_PATTERNS
# echo $PWD

# -L turns on following symbolic links
input=$( find -L $INPUT_FILE_DIR -name "$FIND_PATTERNS" -and -type f -print0 | xargs -0 echo )
read -r -a input_files_array <<< $input

mkdir -p "$OUT_DIR"

# Loop through all the samples
for i in "${!input_files_array[@]}"; do
	# /home/rada/Documents/CGL/JASEN/exp/make_validation_summary/results/Acinetobacter_baumannii_MR130386-1/AMRFinderPlus/aMRFinderPlus.tsv
	FULL_FILE_PATH="$PWD"/$(echo ${input_files_array[$i]})
	# Acinetobacter_baumannii_MR130386-1
	SAMPLE_DATA=$(echo "$FULL_FILE_PATH" | awk -F/ -v sample_info=$SAMPLE_DATA_IN_FULL_PATH_LOC '{print $sample_info}')
	# Parse sample data into species and sample name
	SPECIES=$(echo $SAMPLE_DATA | awk -F_ '{ printf "%s_%s", $1, $2 }')
	SAMPLE_NAME=$(echo $SAMPLE_DATA | awk -F_ '{ printf "%s", $3 }')
	cp "$FULL_FILE_PATH" "$OUT_DIR"/"$SPECIES""_""$SAMPLE_NAME""_""$FILE_TYPE"
done