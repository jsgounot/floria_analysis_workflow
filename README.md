# Glopp analysis workflows

## System requirements

Workflows are based on [Snakemake](https://snakemake.readthedocs.io/en/stable/) and [Conda](https://docs.conda.io/en/latest/) under a Linux environment:

* Ubuntu 18.04.1 LTS
* Conda 4.11.0
* Snakemake 6.15.5

Most of the softwares used in the pipeline will be automatically downloaded by snakemake with conda during first launch following recipies found in `workflow/envs`. Only exceptions are Strainberry and Glopp (**need to added to bioconda once we're done**) which need to be hardcoded in `workflow/strainberry.snk` and  `workflow/glopp.snk`.

Pipeline should take around ~ 24h with 16 CPUs. Runtime is not linear to CPUs count since some part of the pipeline such are variant callers are not multi-threaded.

## Synthetic dataset

The synthetic dataset contains 3 experiments which can be found in `synthetic.json`. With provided config file, sequences will be automatically downloaded from NCBI.

Recommanded command line to run the pipeline:

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

One last example could be `operams.strainberry.nanopore`, which means strainberry with nanopore reads, based on a mapping with operams assemblies (and not the mapping reference from the config file).

## Hardcoded options

Because all potential combinations cannot be informed in file paths, some options were hardcoded for this workflow but might not fit perfectly with other data. For example, the estimated Glopp error rate is uniquely defined for each mode, you will need to modify the snakemake file to change these values.
