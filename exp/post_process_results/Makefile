SHELL = /bin/bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

all:
	python3 gather_files.py -p '**/mlst.json' -s '.json' /data/CGL/JASEN/results /data/CGL/JASEN/mlst_json/