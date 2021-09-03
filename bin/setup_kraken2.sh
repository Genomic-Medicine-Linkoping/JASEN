#!/bin/bash

# http://linuxcommand.org/lc3_man_pages/seth.html
#set -x
# http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
#IFS=$'\n\t'

cd /data/CGL/JASEN
kraken2-build --download-library bacteria --db k2db --threads 90 --use-ftp
kraken2-build --download-taxonomy --db k2db --threads 90 --use-ftp
cd k2db/taxonomy
# These go to /data/CGL/JASEN/k2db/taxonomy directory instead of the ones 
# attempted to be downloaded with kraken2-build --download-taxonomy --db k2db ...
wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz
kraken2-build --build --db k2db --threads 1
