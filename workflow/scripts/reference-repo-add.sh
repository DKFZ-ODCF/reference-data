#!/usr/bin/env bash

# to be called by Snakemake with the 'onsuccess'-trigger to create a git repository from the created repository
set -euf -o pipefail

echo $1 $2 $3
DATA_ROOT_DIR=$1
SPECIES_DIR=$2
ASSEMBLYDIR_VARIANT=$3

snakemake -s ../../refdata-split/workflow/Snakefile --use-conda --configfile config.json --config genome_config=saccharomyces_cerevisiae.json ensembl_release=98 --jobs 2 --latency-wait 20 \
  --archive "${DATA_ROOT_DIR}/${SPECIES_DIR}/${ASSEMBLYDIR_VARIANT}/workflow.tar.gz" #FIXME workflow files not included in archive

snakemake -s ../../refdata-split/workflow/Snakefile --use-conda --configfile config.json --config genome_config=saccharomyces_cerevisiae.json ensembl_release=98 --jobs 2 --latency-wait 20 \
--report "${DATA_ROOT_DIR}/${SPECIES_DIR}/${ASSEMBLYDIR_VARIANT}/report.html"

cd ${DATA_ROOT_DIR} &&
echo "now in $(pwd)" &&
#cd ${SPECIES_DIR} &&
git status && repo_exists=1 || git init &&
git lfs track "*.tar.gz"
git checkout -b "${SPECIES_DIR}/${ASSEMBLYDIR_VARIANT}-dev"

git add "${SPECIES_DIR}/${ASSEMBLYDIR_VARIANT}/workflow.tar.gz"
git commit -m 'add workflow archive'
git add "${SPECIES_DIR}/${ASSEMBLYDIR_VARIANT}/report.html"
git commit -m 'add report'