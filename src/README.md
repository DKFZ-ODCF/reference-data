# Common data preparation

## Using Snakemake 

Idea how to use Snakemake for preparing data in the reference-data repository. 

Write a versatile `Snakefile` with corresponding wrappers for common tasks. 

The data for each genome is stored in `JSON`-files, that could be retrieved via the ENSEMBL API. 

## Configuration

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
   
### How to run

1. If Snakemake is not installed, create environment
```
conda env create -f src/environment.yml
conda activate reference-data
```
2. Adjust `basedir` in `src/config.json` to installation directory. 
 - [ ] Make seperate directories for code (`src`) and data. 
 
3. Obtain copy of Snakemake-wrappers and adjust `wrapperdir` in `src/config.json`
```
git clone https://github.com/mobilegenome/snakemake-wrappers.git
```
4. Request genome-information from ENSEMBL and execute Snakemake
```bash
assembly_name="heterocephalus_glaber_female"
curl "https://rest.ensembl.org/info/genomes/${assembly_name}?" -H "Content-type:application/json" > cache/"${assembly_name}.json" 
snakemake -s src/Snakefile --configfile cache/"${assembly_name}.json"  
```

 

### TODOs
- Snakemake allows defining remote-input (see https://snakemake.readthedocs.io/en/stable/snakefiles/remote_files.html) and incorporation into the dependency-graph. 
However, I have not found a way to configure this for the DKFZ proxy server. 

- Additional Snakemake-rules will then process the sequence files further and create additional files, indeces, etc.
  
