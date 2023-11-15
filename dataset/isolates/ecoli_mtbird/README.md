#### Download ecoli reads and assemblies

This snakemake will download Illumina and Nanopore reads, as well as genbank assemblies linked to E.coli isolates used in Floria paper.

1. Install the conda environment `ncbitools.yaml` if not done already

```
mamba env create -n ncbitools --file ncbitools.yaml
```

2. Use the snakemake pipeline

```
snakemake -s download.snk --use-conda -c {cores} -p --resources ncbi_limit=1
```