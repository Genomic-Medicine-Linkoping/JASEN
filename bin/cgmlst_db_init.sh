#!/bin/bash
set -e
set -u
set -o pipefail

# https://stackoverflow.com/a/16496491
usage() { 
    echo -e "Usage: $0
[-c <NUM_CPUS>]
[-d <PATH_TO_WHERE_DB_SHOULD_BE>]
[-u <CGMLST_SCHEMA_URL>]" 1>&2
    exit 1
}

while getopts ":c:u:d:" o; do
    case "${o}" in
        c)
            c=${OPTARG}
            #((s == "Acinetobacter_baumannii" || s == "Enterobacter_cloacae" || s == "Enterococcus_faecalis" || s == "Enterococcus_faecium" || s == "Escherichia_coli" || s == "Klebsiella_pneumoniae" || s == "Mycobacterium_abcessus" || s == "Mycobacterium_avium" || s == "Mycobacterium_intracellulare" || s == "Mycobacterium_malmoense" || s == "Mycobacterium_tuberculosis" || s == "Pseudomonas_aeruginosa" || s == "Staphylococcus_argenteus" || s == "Staphylococcus_aureus" || s == "Stenotrophomonas_maltophilia" || s == "Streptococcus_pyogenes" || s == "Citrobacter_braakii" || s == "Clostridioides_difficile" || s == "Corynebacterium_striatum" || s == "Enterococcus_gallinarum" || s == "Klebsiella_aerogenes" || s == "Mycobacterium_chimera" || s == "Mycobacterium_malmoense" || s == "Mycobacterium_kansasii" || s == "Mycobacterium_chelonae" || s == "Mycobacterium_celatum" || s == "Mycobacterium_marinum" || s == "Mycobacterium_szulgai" || s == "Mycobacterium_gordonae" || s == "Mycobacterium_scrofulaceum" || s == "Mycobacterium_xenopi" || s == "Mycobacterium_africanum" || s == "Mycobacterium_bovis" || s == "Proteus_mirabilis" || s == "Proteus_vulgaris" || s == "Salmonella_enterica")) || usage
            ;;
        u)
            u=${OPTARG}
            ;;
        d)
            d=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${c}" ] || [ -z "${u}" ] || [ -z "${d}" ]; then
    usage
fi

mkdir -p ${d}
wget -c ${u} -O chewiedb.zip
unzip -o chewiedb.zip -d ${d} 
chewBBACA.py PrepExternalSchema -i ${d} -o ${d}/schema --cpu ${c}
touch database.rdy

