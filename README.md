# Glopp analysis workflows

## System requirements

Workflows are based on [Snakemake](https://snakemake.readthedocs.io/en/stable/) and [Conda](https://docs.conda.io/en/latest/) under a Linux environment:

* Ubuntu 18.04.1 LTS
* Snakemake 7.18.1

Pipeline should take around ~ 24h with 16 CPUs. Runtime is not linear to CPUs count since some part of the pipeline such are variant callers are not multi-threaded. 

## Conda environments and software paths

Your main environment is supposed to have biopython and pigz installed as well. On ubuntu:

```bash
conda install -c bioconda biopython
sudo apt-get install pigz
```

Most of the other softwares used in the pipeline will be automatically downloaded by snakemake with conda during first launch, following recipes found in `workflow/envs`. If you want to use a different conda environments, you can replace associated environment names in the `condaenvs.json` file.

For people with files number quota, a merged version of most of the environments `merged.yaml` (excluding `operams`, `strainberry` and `strainxpress`) is available in `workflow/envs/` and can be used to replace most of the environments yaml files in the `condaenvs.json`.

**IMPORTANT**: For some softwares, you will need to install them locally: [strainberry](https://github.com/rvicedomini/strainberry), [strainxpress](https://github.com/kangxiongbin/StrainXpress), [glopp](https://github.com/bluenote-1577/glopp) and [opera-ms](https://github.com/CSB5/OPERA-MS). You can install just the ones you need, and add the executable path in `softpaths.json`.

## Testing the pipeline

You can test the pipeline with this command line.

```bash
snakemake -s single_species_synthetic.snk --configfile ctest.json --use-conda --conda-prefix ./conda --cores 4 -d ./res
```

This will run the pipeline with a test of two E. coli samples at low coverage.

## Synthetic dataset

The synthetic dataset contains 3 experiments which can be found in `synthetic.json`. With provided config file, sequences will be automatically downloaded from NCBI.

Recommended command line to run the pipeline:

```bash
snakemake -s single_species_synthetic.snk --configfile synthetic.json --use-conda --cores 24 --resources ncbi_load=1 --attempt 3
```

`--resources ncbi_load=1` option is a safegard to limit the number of NCBI call to 1 and <u>should not be removed</u>. 

**Adding new sequences**: You have three options associated with a json key, at least one of them should be provided for each sample and keys are evaluated with this order:

* <u>refpath</u>: a path to a local file.
* <u>ncbiasbly</u>: a NCBI assembly identifier. Note that one assembly can contains multiple sequences (contigs, plasmids) which will be concatenated in the process.
* <u>ncbinuc</u>: A single NCBI identifier like a RefSeq identifier. Ideal if you just want one contig.

## Understanding results file path naming

One phasing is the combination of multiple software with different options. In snakemake, this variations are directly hardcoded in the result directory path, as a requiered input for the `all` rule of the main snakemake script. 

For example, `stats/assemblies/ecoli/glopp.inpref.hybrid.longshot.nanoprep1/mummer/circos/done.empty` contains `glopp.inpref.hybrid.longshot.nanoprep1`. This indicates that we used glopp as phasing algorithm, with input reference (`inpref`), a hybrid phasing based on longshot mapping and a post-phasing assembly method called nanoprep1.

Another example could be: `strainxpress.fast` which simply mean the phasing is done by strainxpress with the fast option. Strainxpress does not use mapping and variant calling to run.

One last example could be `operams.strainberry.nanopore`, which means strainberry with nanopore reads mapping against OPERA-MS assembly.

## Hardcoded options

Because all potential combinations cannot be informed in file paths, some options were hardcoded for this workflow but might not fit perfectly with other data. For example, the estimated Glopp error rate is uniquely defined for each mode, you will need to modify the snakemake file to change these values.

## Known issues

* Wtdbg2 and Strainberry (that uses Wtdbg2) crash 'randomly': This happens on some configuration and sadly there is nothing to do about it. The only way is to use the `retries` option (usually set to 3 on my tests).
* StrainXpress uses too much memory. One way is to limit the number of CPU when calling strainxpress.