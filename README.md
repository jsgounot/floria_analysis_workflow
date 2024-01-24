## Floria pipelines

This repository contains 1) all workflows used to run and compare phasing solutions, including regular assembly software, split phasing and real reads approaches, 2) production workflow to directly run and process you metagenomic reads, and 3) all datasets used in the inital Floria paper.

### Conda environments and software paths

#### System

Workflows have been written using [Snakemake](https://snakemake.readthedocs.io/en/stable/) and [Conda](https://docs.conda.io/en/latest/) under a Linux environment:

* Ubuntu (18.04.1 LTS)
* Snakemake (7.18.1 and 7.32.4)

#### Conda

Your main conda environment is supposed to have snakemake, biopython and pigz installed as well. On ubuntu:

```bash
mamba install -c conda-forge -c bioconda biopython snakemake
sudo apt-get install pigz
```

Most of the other softwares used in the pipeline will be automatically downloaded by snakemake with conda during the first launch, following recipes found in `workflow/envs`. If you want to use different conda environments, you can replace associated environment names in the `condaenvs.json` file.

#### Softwares

Some softwares are not available on conda and need to be installed locally: [strainberry](https://github.com/rvicedomini/strainberry), [strainxpress](https://github.com/kangxiongbin/StrainXpress), [floria](https://github.com/bluenote-1577/glopp), and [opera-ms](https://github.com/CSB5/OPERA-MS). You can install just the ones you need, and add the executable paths in `softpaths.json`.

For Strainberry, I recommend installing my small fork version of the software [here](https://github.com/jsgounot/strainberry) that is more suitable for snakemake pipeline (see details in [this commit](https://github.com/rvicedomini/strainberry/commit/153a84cedb2ed07590af5a6ba0e01899389de1eb)).

If you plan to use the kraken classification approach to identify reference sequences, you'll need to specify both the kraken database (`krakendb`) and the **gzipped** kraken database library fasta file (`krakendb_lib`). If not specified, the pipeline will assume the library file is found within `path/to/krakendb/library/library.fna.gz`. For gut microbiome, you can use the [UHGG database](http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.1/kraken2_db_uhgg_v2.0.1/). Don't forget to compress the fna file within the library directory.

### Production pipeline

If you want to process your metagenomic samples, please use `multiple_species_production.snk`.

Here is an example of the configuration file you can use:

```json
{
    "outputs": [
        {
            "assembly": {
                "name": "kraken_ref",
                "mode": "nanopore"
            },
            "group": "all",
            "phasing": {
                "assembler_preset": "none",
                "assembler_rtype": "long_reads",
                "fmode": "none",
                "name": "floria",
                "post_assembler": "wtdbg2",
                "readtype": "nanopore",
                "assembler_preset": "nano"
            },
            "vcalling": "longshot"
        }
    ],
    "samples": {
        "Sample_name_1": {
            "nanopore": "path/to/your/sample/reads.fastq.gz",
        }
    }
}
```

Please check [the related section in the configuration files guide](config_desc.md#variant-calling) for more information.

### Simulation and assessment pipelines

This repository contains multiple main pipelines (or workflows), each of them being a snakemake file `.snk` in the root folder and composed of multiple subworkflows you can find in the `workflow` folder. 

#### Single species synthetic

`single_species_synthetic.snk`

Produce synthetic reads of multiple strains **from the same species** to produce phasing. For mapping-based approach, reads are mapped against one user-designated reference genome. This is the simplest way of testing phasing algorithms. Note that you can have multiple species that will be processed independently with the same run.

#### Single species subsampling

`single_species_subsampling.snk`

Subsample real reads of multiple strains **from the same species** to produce phasing. Excluding this initial part, the pipeline is similar to *single species synthetic* one.

#### Multiple species / Metagenome synthetic

`multiple_species_synthetic.snk`

Produce synthetic reads **from multiple strains and species to simulate a metagenome** and phase without prior-knowledge of existing species. For reference-based approaches such as Floria or Strainberry, kraken-based or merged assembly are available (see below).

### Note on the reference base approaches for multiple species (including production) pipelines

For multiple species, two approaches are available to generate reference sequence(s) that will be used by Floria or Strainberry.

**Option A**: All reads are assembled within on single reference genome that can be used as reference. This would work with Floria but will lead to very bad results for Strainberry, that assume one single species.

**Option B**: Reference species are defined using a Kraken based approach where reads are classified against a Kraken database and only species with estimated coverage higher than *n*X (default = 5X) are kept. Reads are then assigned against a unique species, based on similarity threshold; and each species is processed individually with Floria or Strainberry. More details are available within our paper.

### Configuration file

A configuration file is requiered to run the pipeline, and examples of configuration files can be found within the `config` folder. A description of configuration file is defined in [a dedicated markdown file](config_desc.md).

### Running the pipeline

The pipeline can be run like this (dry-mode):

```bash
snakemake -s {pipeline_name}.snk --configfile {config.json} -d {outdir} --use-conda --cores {cores} -n
```

I recommend some other options:

* `--conda-prefix {absolute_path}` to root your conda environment within a specific folder.
* `--attempt 2` to be sure that random errors (for example memory) does not make your whole pipeline crash
* ` --resources ncbi_load=1` to ensure that you limit requests against NCBI server (**important**)
* `--rerun-incomplete` to rerun rules that might have failed before
* `--rerun-triggers mtime` to avoid to unnecessary rerun some samples
* `--scheduler greedy` do not use the usual DAG scheduler, can make your pipeline to execute faster
* `--keep-going` to avoid that the pipeline stops at the first error

Some options might be useful in specific conditions

* `-n` to dry-run 
* `--notemp` if you want to keep all temporary files (including large bam files)
* `-p` to show used command lines

#### CPU usage

If benchmarking is not a concern for you, I also recommend to slightly specify more cores than what you configuration offer, fore example 50 CPUs instead of 46. Advantages are multiples: Some used softwares are only partially using all their CPU. Some rules requiere very little CPU usage but can still take some time to be completed, this avoid to have unecessary bottlenecks. Note that increasing this number too much can also lead to some weird and untracktable errors (WTDBG2 for example), so please be cautious. 

### Outputs

#### Production pipelines

The main ouputs are in the `results` folder, including for Floria a ploidy file (`tsv.gz`) and an assembly file containing the concatenation of all haplosets (`fa.gz`).

Ploidy files contain all information related to each phased contig, first column are copied from Floria output ([description here](https://phase-doc.readthedocs.io/en/latest/usage.html)), extra columns are defined here:

| Field               | Description                                                  |
| ------------------- | ------------------------------------------------------------ |
| contig_size         | Contig size                                                  |
| taxid               | Taxonomic ID from the Kraken database                        |
| kraken_abundance    | Kraken abundance for this taxid                              |
| kraken_#covered     | Number of reads covering the taxid                           |
| kraken_#assigned    | Number of reads assigned to this taxid                       |
| kraken_rank         | Kraken taxid rank                                            |
| kraken_name         | Kraken taxid name                                            |
| kraken_assigned     | Should be similar to `kraken_#assigned`, for debug           |
| kraken_readlength   | The cumulative reads length of assigned reads                |
| kraken_assigned_abu | Should be similar to `kraken_abundance`, for debug           |
| kraken_reflen       | Kraken DB reference size                                     |
| kraken_cov          | Coverage estimation based on reads cumulative length and ref size |
| kraken_reflen_mb    | Kraken ref size in Mb                                        |
| bdepth_npos         | Number of position from the extracted reference used for mapping, for debug |
| bdepth_sum_depth    | Total cumulative depth of all position after mapping         |
| bdepth_mean_cov     | Mean coverage of the reference                               |
| af_nba              | Number of SNPs after variant calling                         |
| af_con              | Genome wide average of allele frequencies with confident alleles = `AC` / (`DP` - `AM`) |
| af_amb              | Genome wide average of allele frequencies with ambiguous alleles =  `AC` / `DP` |

An estimation of the strains count <u>at the species level</u> can be obtained from this file:

```python
import numpy as np
import pandas as pd

fname = 'path/to/your/file.tsv.gz'
df = pd.read_csv(fname, sep='\t')

# Minimum 15HAPQ used here!
use_col = 'average_straincount_min15hapq'

# Defining a minimum strain count of 1
df[use_col + '_oneplus'] = np.where(df[use_col] > 1, df[use_col], 1)

# Strain count normalized by contig size
df['alpw'] = df[use_col + '_oneplus'] * df['contig_size']
columns = ['alpw', 'contig_size', 'bdepth_npos', 'bdepth_sum_depth']
mdf = df.groupby('taxid')[columns].sum().reset_index()

# Adding coverage information
mdf['mean_cov'] = mdf['bdepth_sum_depth'] / mdf['bdepth_npos']
mdf['alp'] = mdf['alpw'] / mdf['contig_size']

# Adding Kraken information
kabu = df[['taxid', 'kraken_abundance']].drop_duplicates()
mdf = mdf.merge(kabu, on='taxid', how='left')

# Rounded value = Approximation of the number of strains you have here
mdf['alp_r'] = mdf['alp'].round()

mdf.head()
```

Please check paper methodology if you want more details regarding this part.

#### Assessment pipelines

Results for each main *phases* of the pipeline are saved into different folders. The main outputs are:

* `assembly` contains the link to all prior assemblies (such as flye or megahit)
* `phasing` contains all the phasing results, including the concatenated haplotypes 
* `stats` contains haplotype statistics (for assessement pipelines) and circos files
* `refcomp/*.fastani.txt` & `refcomp/*.mummer.txt` for the reference comparisons when enabled
* `benchmarks` contains the individual benchmark of the most CPU intensive rules

#### Notes on benchmarking

While most processes use CPU quite efficiently, some software underuse the number of available CPUs (e.g strainberry during the longshot process, using only one CPU). This results in high CPU wall time but low CPU loading, both being provided by snakemake (respectively `cpu_time` and `mean_load`). While it is tempting to only report CPU loading, we find this metric unfair for CPU efficient softwares. For Floria's paper, we decided to use snakemake CPU wall time normalized by the number of CPU used for the specific rule. While this can't be done directly by snakemake and requires the user to remember the number of CPUs used for each rule, this value provides a more accurate representation of each rule's CPU efficiency.

### Citations

TBD
