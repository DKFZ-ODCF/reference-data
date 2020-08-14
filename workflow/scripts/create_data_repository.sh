#!/usr/bin/env bash

# to be called by Snakemake with the 'onsuccess'-trigger to create a git repository from the created repository
set -euf -o pipefail

# Test if git and git annex are installed
git --version  && echo "git installed"  || (echo "git not installed"; exit) &&
git-annex version && echo "git-annex installed" || (echo "git annex not installed"; exit)

ASSEMBLYDIR_VARIANT=$1

cd ${ASSEMBLYDIR_VARIANT} &&
conda env create -f workflow/envs/git-annex.yml &&
conda activate git-annex &&
git init &&
git annex init &&
git annex add . &&
git commit -m 'add files'