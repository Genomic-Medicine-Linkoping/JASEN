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

while getopts ":s:p:" o; do
    case "${o}" in
        s)
            s=${OPTARG}
            ((s == 45 || s == 90)) || usage
            ;;
        p)
            p=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${s}" ] || [ -z "${p}" ]; then
    usage
fi

mkdir -p ${d}
wget -c ${u} -O chewiedb.zip
unzip -o chewiedb.zip -d ${d} 
chewBBACA.py PrepExternalSchema -i ${d} -o ${d}/schema --cpu ${c}
touch database.rdy

