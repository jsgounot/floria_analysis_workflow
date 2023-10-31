## Floria pipelines

This repository contains all workflow used to run and compare phasing solutions, including regular assembly software, split phasing and real reads approaches.

### Conda environments and software paths

#### System

Workflows have been developed using [Snakemake](https://snakemake.readthedocs.io/en/stable/) and [Conda](https://docs.conda.io/en/latest/) under a Linux environment:

* Ubuntu 18.04.1 LTS
* Snakemake 7.18.1

#### Conda

Your main environment is supposed to have biopython and pigz installed as well. On ubuntu:

```bash
conda install -c bioconda biopython
sudo apt-get install pigz
```

Most of the other software used in the pipeline will be automatically downloaded by snakemake with conda during the first launch, following recipes found in `workflow/envs`. If you want to use different conda environments, you can replace associated environment names in the `condaenvs.json` file.

A merged version of most of the environments `merged.yaml` (excluding `operams`, `strainberry` and `strainxpress`) is available in `workflow/envs/` and can be used to replace most of the environments yaml files in the `condaenvs.json`.

#### Softwares

You need to install some software locally: [strainberry](https://github.com/rvicedomini/strainberry), [strainxpress](https://github.com/kangxiongbin/StrainXpress), [floria](https://github.com/bluenote-1577/glopp), and [opera-ms](https://github.com/CSB5/OPERA-MS). You can install just the ones you need, and add the executable path in `softpaths.json`.

For Strainberry, I recommend installing my small fork version of the software [here](https://github.com/jsgounot/strainberry) that is more suitable for snakemake pipeline (see details in [this commit](https://github.com/rvicedomini/strainberry/commit/153a84cedb2ed07590af5a6ba0e01899389de1eb)).

### Pipelines description

This repository contains multiple main pipelines (or workflows), each of them being a snakemake file `.snk` in the root folder and composed of multiple subworkflow you can find in the `workflow` folder. 

#### Single species synthetic

Produce synthetic reads of multiple strains **from the same species** to produce phasing. For mapping-based approach, reads are mapped against one user-designated reference genome. This is the simplest way of testing phasing algorithms. Note that you can have multiple species that will be treated independently in the same configuration file.

#### Single species subsampling

Subsample real reads of multiple strains **from the same species** to produce phasing. Excluding this initial part, the pipeline is similar to *single species synthetic* one.

#### Multiple species / Metagenome synthetic

Produce synthetic reads **from multiple strains and species to simulate a metagenome** without knowledge of existing species. For reference-based approaches such as Floria or Strainberry, this means we need to first identify the potential species within the metagenome. This is done using a kraken-based approach.

#### Multiple species spike-in

A slightly different and minor approach of the *multiple species synthetic* pipeline where simulated reads are spike-in within a real metagenomic reads dataset. Mainly used to assess the ability of the split-kraken part to correctly identify species and phasing algorithm to reconstruct strains within such samples.

#### Multiple species production (Real reads)

Similar to *multiple species synthetic* but for real metagenomic samples, without the synthetic reads part. This pipeline can be used for your sample but most likely will require tuning some parameters first. 

### Testing the pipeline

You can test the pipeline with this command line.

```bash
snakemake -s single_species_synthetic.snk --configfile configs/single_species_synthetic_test.json -d ./res/single_species_synthetic_test --use-conda --conda-prefix ./conda --cores 16 --resources ncbi_load=1 --attempt 3
```

This will run the pipeline with multiple strains from *E. coli* and *K. pneumonia*. 

### Synthetic dataset generation

#### Reference genomes

You can either provide a local path to a reference genome or download the reference on NCBI, using the following keywords in the configuration file:

* <u>refpath</u>: a path to a local file.
* <u>ncbiasbly</u>: a NCBI assembly identifier. Note that one assembly can contain multiple sequences (contigs, plasmids) which will be concatenated in the process.
* <u>ncbinuc</u>: A single NCBI identifier like a RefSeq identifier. Ideal if you just want one contig.

See the `configs/single_species_synthetic_test.json` for an example. Snakemake option `--resources ncbi_load=1` is a safeguard to limit the number of NCBI call to 1 and <u>should not be removed</u> if you download NCBI genomes.

#### Reads simulation

Illumina reads are automatically simulated with [ART](https://github.com/scchess/Art), nanopore reads with [badread](https://github.com/rrwick/Badread), pacbio reads with [pbsim](https://github.com/yukiteruono/pbsim3). For nanopore reads, the configuration file options `circular` (default: `true`) and `nanopore` (default: `average`) can be used to define circularity and nanopore read quality. For all read simulators, a seed can be defined to ensure reproducibility. 

### Path name presets and combinations

Snakemake is fundamentally relying on path name and [wildcards](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards) to define how to produce an output file. Since a phased assembly can result from multiple chained operations, they must be reflected within the final path.  

The path base of an output file is `stats/assemblies/{group}/{basename}/mummer/circos/done.empty` where `group` refers to one of the main groups in your configuration file and `basename` encodes how to make the phasing. A very simple `basename` can be any assembler or reference-free phasing software such as `megahit` or `strainxpress.regular` (regular mode for StrainXPress). 

For reference based approaches, `basename` *can* be composed of multiple fields based on this order: `{assembling}.{ass_reads}.{vcalling}.{phaser}.{phaser_mode}.{subassembler}.{subassembler_reads}.{subassembler_preset}`. 

#### Fields description

##### Assembling

* `inpref`: For *single species* approach, will take the reference genome indicated in the configuration file.
* `{assembler_name}`: For *single species* approach, will take the assembly generated by `{assembler_name}`
* `kraken_ref`: For *multiple species* approach, will take a reference genome and bam file generated using the kraken-mapping approach with each species processed individually. Requiere `{ass_reads}` that indicates which reads will be used for the kraken analysis.
* `kraken_presplit`: Similar to kraken_ref, but all species are not split but processed in the same batch.

##### Vcalling

Informs which variant caller to use, can be either `longshot`, `lofreq` or a custom handmade variant caller called `binomial_custom`. We do not recommend using the latest as it oversimplifies the variant calling process. `longshot` should be the default choice, `lofreq` does ***not*** work well with nanopore long reads.

##### Phaser and phaser mode

Can be either `floria` or `strainberry` for now.

Phaser mode is only used for *single species* approach with floria and defines whether you want to use `illumina`, `nanopore`, `pacbio` or an `hybrid` of illumina and nanopore. Note that for `pacbio` you also need to define the preset for read simulation.

##### Sub assembler

For floria only, define if you want to use `abysspe` or `wtdbg2`, using either (`{subassembler_reads}`) the `short_reads` or the `long_reads` of glopp output.  For `wgtdbg2` you also need to add the `subassembler_preset` corresponding for now to `nano` (`-x preset2 -e 5 -l 1000 -L 3000 -S 1 -R'`) for nanopore and `hifi` (`-x ccs -R'`) for hifi reads.

### Citations

TBD
