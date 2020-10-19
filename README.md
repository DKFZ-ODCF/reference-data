# Reference data 

![CI](https://github.com/DKFZ-ODCF/reference-data/workflows/CI/badge.svg?branch=development&event=push)
[![Snakemake](https://img.shields.io/badge/snakemake-≥5.20-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)


A snakemake-based pipeline to set up reference genome and corresponding data.

## Concept

Idea how to use Snakemake for preparing data in the reference-data repository. 

Write a versatile `Snakefile` with corresponding wrappers for common tasks. 

The data for each genome is stored in `JSON`-files, that could be retrieved via the ENSEMBL API. 

### Data management

Datasets created by this pipeline will the follow a four level directory structure with the pattern

`$species/$primary_assembly_variant/$sequence_variant/$snv_variant`

Here, `$species` refers to the scientific name of the organism,
 `$primary_assembly` to the assembly variant provided by public resources, 
 `$naming_variant` to sequence differences created by addition / removal of auxillary sequence (e.g. *PhiX*), and
 `$snv_variant` to variants caused by replacing nucleotides. 
 
For brevity, the levels are also referred to as 1-4.

### Versioning data

Each dataset will be initialized as an own `git` repository so that workflows 
and analyses can always refer to versions of reference data used.

### Using 'git annex'

Reference data is initialized as a 'git' repository with all files added via `git annex`. 

As a user, you can retrieve the data as following: 
```
git clone <data-repository> <target-dir>
cd <target-dir>    
git annex get # retrieve file
```

## Configuring the pipeline

A global configuration is stored in `src/config.json`.
This defines the base-directory of the repository and the ENSEMBL version to be used for downloading.
*Note* the ENSEMBL release version used may vary between species/genomes depending on the researches demands. For all annotation files, the release version is part of the filename. 

## Routines

 - [x] Download genome sequence, including CHECKSUM verification 
 - [x] Download annotation, including CHECKSUM verification 
 - [x] Un(g)zip data
 - [x] bwa index 
 - [x] STAR index (not tested)
 - [ ] add sequences (e.g. PhiX)        
 - [ ] replace chromosome prefixes
   
## How to run

1. If Snakemake is not installed, create environment
```
conda env create -f workflows/envs/snakemake-base.yml
conda activate snakemake-base
```
2. Configuration
  - Adjust the directories in `config/config.json` .  
 
 -  The pipelines uses a forked version of the snakemake-wrapper repository. In order to ensure reproducibility, change the
     `wrapperdir" and `wrapper_custom_commit` variables to the commit you will be using. For example: 
     ```
    "wrapper_dir": "https://raw.githubusercontent.com/DKFZ-ODCF/snakemake-wrappers",
    "wrapper_custom_commit": "6b32a15a"
    ```
     
4. Request genome-information from ENSEMBL and execute Snakemake

```bash
    assembly_name="heterocephalus_glaber_female"
    curl "https://rest.ensembl.org/info/genomes/${assembly_name}?" -H "Content-type:application/json" > cache/"${assembly_name}.json" 
    snakemake --use-conda --configfile resources/examples/drosophila_melanogaster.json --config configfile=config/config.json ensembl_release=98 --jobs 4 --latency-wait 20
```


## Snakemake wrappers

Whenever possible, the workflow uses wrapper to externalize tasks. In order to allow some independence from the offical *snakemake-wrappers* repository (LINK)
that provides universally applicable wrappers, we're using a [forked wrappers repository](https://github.com/DKFZ-ODCF/snakemake-wrappers)

## Development

### Testing

This repository is using github *actions* for automated testing. The workflow described in `workflows/test.yaml` performs 
a dry run of the pipeline, a run to create a summary report for inspection of the snakemake workflow and a actual test run using 
the *Drosophila melanogaster* genome. 

The test workflow is triggered at each push of the development branch.  
