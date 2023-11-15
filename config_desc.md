## Description of the configuration files

A configuration file is requiered to run the pipeline. The configuration file is split into 3 different sections:

1. The `samples`, information relative to your samples, e.g reads path.
2. The `outputs`, which assembler or phasers to use, which options.
3. The `miscs`, including default number of CPUs or presets.

Since two runs can be very different from one experiment to another, the configuration file cover a lot of usages with **some elements are requiered by some pipelines and not by others**. This might be quite difficult to handle at first, please check the examples in the `configs` folder if you prefere more practical examples. You can also check out some configurations for the most basic usages at the end of this file.

The json configuration file must be like:

```json
{
  "samples": {},
  "outputs": [],
  "miscs": {}
}
```

### The samples

Each first key must be a **group** that should have two or more **samples** (each sample being a strain). With the exception of the production pipeline, a reference genome is requiered for each sample to be used for assessment.

```json
{
	"groupname":{
		"sample1": {},
		"sample2": {}
	}
}
```

For each sample, you can use those keys:

| Key             | Description                                                  | Note                                                         |
| --------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| `refpath`       | The reference genome path [1]                                | You can assemble reads otherwise                             |
| `ncbiasbly`     | The NCBI assembly ID [1] [2]                                 | Can contains multiple contigs, that will be concatenated.    |
| `ncbinuc`       | The NCBI nuc ID [1] [2]                                      | Ideal if you just want one contig.                           |
| `mapping_ref`   | Is this sample genome used as reference for the group        | Only one per group. Mandatory to have one set to true.       |
| `nanopore`      | Nanopore reads path if used                                  |                                                              |
| `nanopore_qual` | Nanopore quality for badread. Can be either `bad`, `good` or `average`. | Predefined presets, see the snakemake workflow for more details. |
| `illumina_r1`   | Illumina R1 path if used                                     |                                                              |
| `illumina_r2`   | Illumina R2 path if used                                     |                                                              |
| `seed`          | Random seed                                                  |                                                              |
| `quantity`      | Target coverage to achieve for either subsampling or simulation | Integer                                                      |
| `circular`      | `true` or `false`. For nanopore reads simulation with badread. |                                                              |

[1] Can only be one.

[2] If used, I strongly encourage you to add the snakemake option `--resources ncbi_load=1` as a safeguard to limit the number of NCBI call to a single one.

### The outputs

This is the most complicated part. It is defined in multiple subsections that can be optional or not depending of what you want to achieve.

#### Group

Is this output for all groups defined in the samples section or just one. If just one, write the group name otherwise write *all*.

#### Assembly

How is defined the group reference genome. Can be either pre-defined or constructed directly. Is a subdictionnary.

##### **1. Predefined**

Just one key

| Key    | Value    | Note |
| ------ | -------- | ---- |
| `name` | `inpref` |      |

##### **2. Assembly with Flye**

| Key    | Value                    | Note              |
| ------ | ------------------------ | ----------------- |
| `name` | `flye`                   |                   |
| `read` | `nanopore` or `illumina` |                   |
| `mode` | `raw`                    | Limited right now |

##### **3. Assembly with Megahit**

| Key    | Value     | Note |
| ------ | --------- | ---- |
| `name` | `megahit` |      |

##### **4. Kraken approach**

In this approach we defined what reference genome(s) are within your samples using Kraken against an UHGG database.

| Key    | Value                    | Note |
| ------ | ------------------------ | ---- |
| `name` | `kraken_ref`             |      |
| `mode` | `nanopore` or `illumina` |      |

##### 5. StrainXPress

While doing strains phasing, StrainXPress is defined as an assembly since it does not need any assembly to work.

| Key    | Value               | Note |
| ------ | ------------------- | ---- |
| `name` | `strainxpress`      |      |
| `mode` | `slow` or `regular` |      |

#### Variant calling

Which variant caller to use. Direct value, can be:

* `freebayes` - can be used for short-reads instead of longshot
* `longshot` - the default option for nanopore
* `lofreq`  - does not work well for nanopore reads
* `binocustom` - a custom and way (too) simple variant caller using binomial distribution probability.

#### Phasing

##### 1. Floria

| Key                | Value                                                  | Note |
| ------------------ | ------------------------------------------------------ | ---- |
| `name`             | `floria`                                               |      |
| `readtype`         | `nanopore` or `illumina`                               |      |
| `fmode`            | `none` or a specific ID defined in `miscs` section     |      |
| `post_assembler`   | `inference`, `megahit`, `abysspe`, `flye` or `wtdbg2`  |      |
| `assembler_rtype`  | `short_reads` or `long_reads`                          |      |
| `assembler_preset` | `none` or a specific preset defined in `miscs` section |      |

Inference post-assembler is simply infering SNPs from the VCF directly within the reference genome. While being extremely fast, it does not take into consideration potential structural variants and can lead to multiple errors. We do not recommand this option in general.

**Example**

```
"phasing": {
	"name": "floria",
	"readtype": "illumina",
	"fmode": "none",
	"post_assembler": "megahit",
	"assembler_rtype": "short_reads",
	"assembler_preset": "none"
}
```

##### 2. Strainberry

| Key        | Value               | Note                         |
| ---------- | ------------------- | ---------------------------- |
| `name`     | `strainberry`       |                              |
| `readtype` | `nanopore` or empty | To add the `nanopore` option |

#### References filtering

A references filtering module is available to ignore phasing results from region that are too similar to each others (between reference genomes). This can be useful if you want to compare phasing results from two strains that have sections with less than X% divergence, and that would be too difficult to untangle with regular phasers. In practice, each reference genomes are pairwised compared to each other with MumMer and only regions with a certain length and similarity will be rejected.

| Key          | Value           | Note     |
| ------------ | --------------- | -------- |
| `length`     | `none` or int   | e.g 1000 |
| `similarity` | `none` or float | e.g 99.5 |

**Example**

```json
"ref_filtering": {
    "length": "none",
    "similarity": "none"
}
```

#### References comparison

A module to easily compare your input references with either FastANI or pairwise Mummer. Can be usueful to quickly check strains / species similarity, identify potential large structural variants or regions with very high or low divergence. Can only be `mummer` or `fastani`. Results can be found within `result/dir/refcomp/{groupname}.{method}.txt`.

### Miscs

Can be used to defined specific part of the pipeline, except the `GENERIC_THREADS` options,  you probably don't want to change it. 

```
{
	"GENERIC_THREADS": 32,
	"GENERIC_THREADS_SPLIT": 8,
	"FLORIA": {
		"PRESET_NAME": "OPTIONS" 
	},
	"WTDBG2": {
		"PRESET_NAME": "OPTIONS"
	}
}
```

### Examples

Generate phasing assembly with either Floria + WTDBG2, Floria + Flye or Strainberry for two communities, composed to either 3 or 2 *E. coli* strains that will be downsampled to 15 or 25X.

```json
{
    "miscs": {
        "GENERIC_THREADS": 32,
        "GENERIC_THREADS_SPLIT": 8
    },
    "outputs": [
        {
            "assembly": {
                "name": "inpref"
            },
            "group": "all",
            "phasing": {
                "assembler_preset": "none",
                "assembler_rtype": "long_reads",
                "fmode": "none",
                "name": "floria",
                "post_assembler": "flye",
                "readtype": "nanopore"
            },
            "ref_filtering": {
                "length": "none",
                "similarity": "none"
            },
            "vcalling": "longshot"
        },
        {
            "assembly": {
                "name": "inpref"
            },
            "group": "all",
            "phasing": {
                "assembler_preset": "nano",
                "assembler_rtype": "long_reads",
                "fmode": "none",
                "name": "floria",
                "post_assembler": "wtdbg2",
                "readtype": "nanopore"
            },
            "ref_filtering": {
                "length": "none",
                "similarity": "none"
            },
            "vcalling": "longshot"
        },
        {
            "assembly": {
                "name": "inpref"
            },
            "group": "all",
            "phasing": {
                "name": "strainberry",
                "readtype": "nanopore"
            },
            "ref_filtering": {
                "length": "none",
                "similarity": "none"
            },
            "vcalling": "longshot",
          	"refcomp": [
							'fastani', 'mummer'
						]
        }
    ],
    "samples": {
        "ecoli_subgroup_2": {
            "542093": {
                "illumina_r1": "path/to/illumina.542093.r1.fastq.gz",
                "illumina_r2": "path/to/illumina.542093.r2.fastq.gz",
                "mapping_ref": "False",
                "nanopore": "path/to/read.542093.fastq.gz",
                "quantity": 25,
                "refpath": "path/to/ref.542093.fasta.gz",
                "seed": 50
            },
            "786605": {
                "illumina_r1": "path/to/illumina.786605.r1.fastq.gz",
                "illumina_r2": "path/to/illumina.786605.r2.fastq.gz",
                "mapping_ref": "True",
                "nanopore": "path/to/read.786605.fastq.gz",
                "quantity": 25,
                "refpath": "path/to/ref.786605.fasta.gz",
                "seed": 50
            }
        },
        "ecoli_subgroup_3": {
            "542093": {
                "illumina_r1": "path/to/illumina.542093.r1.fastq.gz",
                "illumina_r2": "path/to/illumina.542093.r2.fastq.gz",
                "mapping_ref": "False",
                "nanopore": "path/to/read.542093.fastq.gz",
                "quantity": 25,
                "refpath": "path/to/ref.542093.fasta.gz",
                "seed": 50
            },
            "786605": {
                "illumina_r1": "path/to/illumina.786605.r1.fastq.gz",
                "illumina_r2": "path/to/illumina.786605.r2.fastq.gz",
                "mapping_ref": "True",
                "nanopore": "path/to/read.786605.fastq.gz",
                "quantity": 25,
                "refpath": "path/to/ref.786605.fasta.gz",
                "seed": 50
            },
            "899091": {
                "illumina_r1": "path/to/illumina.899091.r1.fastq.gz",
                "illumina_r2": "path/to/illumina.899091.r2.fastq.gz",
                "mapping_ref": "False",
                "nanopore": "path/to/reads.899091.fastq.gz",
                "quantity": 15,
                "refpath": "path/to/ref.899091.fasta.gz",
                "seed": 50
            }
        }
    }
}
```

### Advanced: Path name presets and combinations

Snakemake is fundamentally relying on path name and [wildcards](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards) to define how to produce an output file. Since a phased assembly can result from multiple chained operations, they must be reflected within the final path. In the previous section, the pipeline will automatically reconstruct the requiered pathname based on the specified requirement. Nevertheless, one can still try to input the specific path directly within the snakemake pipeline within the `all` rule, based on this information:

The path base of an output file is `stats/assemblies/{group}/{basename}/mummer/circos/done.empty` where `group` refers to one of the main groups in your configuration file and `basename` encodes how to make the phasing. A very simple `basename` can be any assembler or reference-free phasing software such as `megahit` or `strainxpress.regular` (regular mode for StrainXPress). 

For reference based approaches, `basename` *can* be composed of multiple fields based on this order: `{assembling}.{ass_reads}.{vcalling}.{phaser}.{phaser_mode}.{subassembler}.{subassembler_reads}.{subassembler_preset}`. 