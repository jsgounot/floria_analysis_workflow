#### Download K. pneumonia and B. lichenformis nanopore reads and produce assembly

Download and assemble nanopore reads used in Floria paper. Use Flye and medaka for assembly.

1. Install the conda environment `ncbitools.yaml` if not done already

You can find the environment within the `ecoli_mtbird` directory.

```
mamba env create -n ncbitools --file ncbitools.yaml
```

2. Install the conda environment for the assembly

If you already have flye and medaka in a conda environment, you can also replace the environment names within the snakemake file. Otherwise you can download those into a separate environment:

```
mamba env create -n isoassembly --file isoassembly.yaml
```

or just add them to an existing environment

```
mamba install -c conda-forge -c bioconda flye medaka
```


3. Use the snakemake pipeline

```
snakemake -s download_assemble.snk --use-conda -c {cores} -p
```