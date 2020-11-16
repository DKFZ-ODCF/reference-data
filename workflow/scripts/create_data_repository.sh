#!/usr/bin/env bash

# to be called by Snakemake with the 'onsuccess'-trigger to create a git repository from the created repository
set -euf -o pipefail

# Test if git and git annex are installed
git --version > /dev/null  && echo "git installed"  || (echo "git not installed"; exit) &&
git-annex version > /dev/null && echo "git-annex installed" || (echo "git annex not installed"; exit)

ASSEMBLYDIR_VARIANT=$1

conda env create -f workflow/envs/git-annex.yml &&
conda activate git-annex &&
cd ${ASSEMBLYDIR_VARIANT} &&
git init &&
#git annex init &&
#git annex add . &&
git commit -m 'add files'

echo "git annex repository created"
echo ""
echo "You may now clone this repository via "
echo "$ git clone"
echo "$ git annex init"
echo "$ git remote add origin $(pwd) (TO BE CHANGED TO REMOTELY ACCESSIBLE)"
echo "and selectively download files via "
echo "% git annex sync origin"
echo "$ git annex get <file>"
echo "============================"
