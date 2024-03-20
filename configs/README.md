### Configuration files

This folder contains configuration files and script to generate configuration files based on data downloaded in the dataset folder.

#### Test configuration files

A test configuration file can be generated for each assessment pipeline with the `make_test_configs.py` script. With exception of `test_single.json` that does not require any input file, you will first need to download and sometimes assemble reads using scripts within the dataset folder. Once done, run the script and try each of those configuration file within your root folder:

`single_species_synthetic.snk`

```
snakemake -s single_species_synthetic.snk --configfile configs/test_single.json --use-conda --conda-prefix {conda_prefix} --cores {cores} -d res/test_single --rerun-incomplete -p --rerun-triggers mtime --scheduler greedy --notemp --keep-going --resources ncbi_load=1 -n 
```

`single_species_subsampling.snk`

```
snakemake -s single_species_subsampling.snk --configfile configs/test_subsampling.json --use-conda --conda-prefix {conda_prefix} --cores {cores} -d res/test_subsampling --rerun-incomplete -p --rerun-triggers mtime --scheduler greedy --notemp --keep-going -n
```

`multiple_species_synthetic.snk`

```
snakemake -s multiple_species_synthetic.snk --configfile configs/test_multiple_synthetic.json --use-conda --conda-prefix {conda_prefix} --cores {cores} -d res/test_multiple_synthetic --rerun-incomplete -p --rerun-triggers mtime --scheduler greedy --notemp --keep-going -n
```

Each test file should run in less than an hour with 48 CPUs.

#### Real case configuration files

**Real reads**

Generat the 3 species (E. coli, B. lichenformis and K. pneumoniae) configuration. **You first need to download the reads and make assemblies** using script in the dataset folder.

```
python real_reads.py
```

You can then try to run the pipeline in the main folder:

```
snakemake -s single_species_subsampling.snk --configfile configs/real_reads.json --use-conda --conda-prefix {conda_prefix} --cores {cores} -d res/real_reads --rerun-incomplete -p --rerun-triggers mtime --scheduler greedy --keep-going
```