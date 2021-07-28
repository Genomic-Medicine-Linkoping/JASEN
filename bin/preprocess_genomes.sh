#!/bin/bash

# http://linuxcommand.org/lc3_man_pages/seth.html
#set -x
# http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
#IFS=$'\n\t'

INPUT="genome_data.tsv"
CONT_NAME=$1
OUT_DIR="assets/ref_genomes"
PT="assets/prodigal_training_files"

# TODO: Try to parallelize this
while IFS= read -r LINE; do
  SPECIES=$(echo "$LINE" | cut -f1 | tr _ " ")
  FN=$(echo "$LINE" | cut -f2)
  URL=$(echo "$LINE" | cut -f3)
  UC=$(basename -s .gz $FN)
  BASENAME=$(basename -s .fna $UC)
  COMP_FILE="$OUT_DIR"/"$FN"
  UNCOMP_FILE="$OUT_DIR"/"$UC"
  TRN="$PT"/"$BASENAME"".trn"

	echo "Downloading $SPECIES genome:"
	echo ""
	wget -vc -O "$COMP_FILE" "$URL"
	echo "Appending md5sum of $SPECIES genome file: $FN to $OUT_DIR/md5sums.txt"
	md5sum "$COMP_FILE" >> "$OUT_DIR"/md5sums.txt
  echo "Uncompressing the genome file: $FN"
  zcat "$COMP_FILE" > "$UNCOMP_FILE"
	echo ""
  echo "Creating prodigal training file: $TRN"
  cat "$UNCOMP_FILE" | singularity exec $CONT_NAME prodigal -p single -t $TRN
done < "$INPUT"