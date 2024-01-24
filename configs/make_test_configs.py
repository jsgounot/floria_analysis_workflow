# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-11-01 15:34:03
# @Last Modified by:   jsgounot
# @Last Modified time: 2024-01-18 15:29:56

import pandas as pd
import glob, json, os

# --------------------------------------------------------
# single_species_synthetic.snk

# snakemake -s single_species_synthetic.snk --configfile configs/test_single.json --use-conda --conda-prefix {conda_prefix} --cores {cores} -d res/test_single --rerun-incomplete -p --rerun-triggers mtime --scheduler greedy --notemp --keep-going --resources ncbi_load=1 -n 

jdata = {
	'samples': {
		"ecoli_seed50": {
			"ecoli_k12": {
				"ncbinuc": "NC_000913.3",
				"mapping_ref": "True",
				"quantity": 20,
				"circular": "true",
				"nanopore_qual": "average",
				"seed": 50
			},
			"ecoli_h5": {
				"ncbinuc": "CP010169.1",
				"quantity": 10,
				"circular": "true",
				"nanopore_qual": "average",
				"seed": 50
			},
			"ecoli_O157": {
				"ncbinuc": "NC_002695.2",
				"quantity": 5,
				"circular": "true",
				"nanopore_qual": "average",
				"seed": 50
			},
			"ecoli_w": {
				"ncbinuc": "NC_017635.1",
				"quantity": 30,
				"circular": "true",
				"nanopore_qual": "average",
				"seed": 50
			}
		}
	},

	'outputs': [
		{
			'group': 'all',
			'assembly': {
				'name': 'inpref'
			},
			'vcalling': 'longshot',
			'phasing': {
				'name': 'floria',
				'readtype': 'nanopore',
				'fmode': 'none',
				'post_assembler': 'wtdbg2',
				'assembler_rtype': 'long_reads',
				'assembler_preset': 'nano'
			},
			'ref_filtering':{
				'length': 'none',
				'similarity': 'none'
			}
		},

		{
			'group': 'all',
			'assembly': {
				'name': 'inpref'
			},
			'vcalling': 'longshot',
			'phasing': {
				'name': 'strainberry',
				'readtype': 'nanopore'
			},
			'ref_filtering':{
				'length': 'none',
				'similarity': 'none'
			}
		}
	],

	'miscs': {
		'GENERIC_THREADS': 8
	}
}

outfile = './test_single.json'
with open(outfile, 'w') as f:
    json.dump(jdata, f, indent=4, sort_keys=True)

# --------------------------------------------------------
# single_species_subsampling

# snakemake -s single_species_subsampling.snk --configfile configs/test_subsampling.json --use-conda --conda-prefix {conda_prefix} --cores {cores} -d res/test_subsampling --rerun-incomplete -p --rerun-triggers mtime --scheduler greedy --notemp --keep-going -n

wd = os.path.dirname(os.getcwd())

fname = os.path.join(wd, 'dataset/isolates/ecoli_mtbird/accessions.tsv')
df = pd.read_csv(fname, sep='\t')
df = df[(df['Chromsome or plasmid'].str.strip() == 'Chromosome') | (df['Chromsome or plasmid'] == 'Chromsome')]
gb = df.set_index('Isolate')['Genbank accession'].str.strip().to_dict()
gb = {str(isolate): gbid for isolate, gbid in gb.items()}

fnames = os.path.join(wd, 'dataset/isolates/ecoli_mtbird/nanopore/*/*.fastq.gz')
fnames = sorted(glob.glob(fnames))
fnames = {os.path.basename(os.path.dirname(fname)): fname for fname in fnames}
assert fnames

subgroups = [
    ('786605', '899091'),
    ('786605', '542093', '899091'),
]

coverages = [25, 20, 15, 30]

jdata = {'samples': {}, 'outputs': {}}
for subgroup in subgroups:
    name = f'ecoli_subgroup_{len(subgroup)}'

    config = {}

    for idx, isolate in enumerate(subgroup):
        refpath = os.path.join(wd, 'dataset/isolates/ecoli_mtbird/genbank',
            isolate, gb[isolate] + '.fasta.gz')

        assert os.path.isfile(refpath)

        read = fnames[isolate]
        assert os.path.isfile(read)

        ir = os.path.join(wd, f'dataset/isolates/ecoli_mtbird/illumina/{isolate}/*.cleaned.fastq.gz')
        ir = sorted(glob.glob(ir))
        assert len(ir) == 2

        config[isolate] = {
            'refpath': refpath,
            'mapping_ref': 'True' if idx == 0 else 'False',
            'nanopore': read,
            "illumina_r1": ir[0],
            "illumina_r2": ir[1],
            'seed': 50,
            'quantity': coverages[idx]
            }

    jdata['samples'][name] = config


jdata['outputs'] = [
	{
		'group': 'all',
		'assembly': {
			'name': 'inpref'
		},
		'vcalling': 'longshot',
		'phasing': {
			'name': 'floria',
			'readtype': 'nanopore',
			'fmode': 'none',
			'post_assembler': 'flye',
			'assembler_rtype': 'long_reads',
			'assembler_preset': 'none'
		},
		'ref_filtering':{
			'length': 'none',
			'similarity': 'none'
		}
	},

	{
		'group': 'all',
		'assembly': {
			'name': 'inpref'
		},
		'vcalling': 'longshot',
		'phasing': {
			'name': 'floria',
			'readtype': 'nanopore',
			'fmode': 'none',
			'post_assembler': 'wtdbg2',
			'assembler_rtype': 'long_reads',
			'assembler_preset': 'nano'
		},
		'ref_filtering':{
			'length': 'none',
			'similarity': 'none'
		}
	},

	{
		'group': 'all',
		'assembly': {
			'name': 'inpref'
		},
		'vcalling': 'longshot',
		'phasing': {
			'name': 'strainberry',
			'readtype': 'nanopore'
		},
		'ref_filtering':{
			'length': 'none',
			'similarity': 'none'
		}
	}
]

jdata['miscs'] = {
	'GENERIC_THREADS': 32,
	'GENERIC_THREADS_SPLIT': 8
}

outfile = './test_subsampling.json'
with open(outfile, 'w') as f:
    json.dump(jdata, f, indent=4, sort_keys=True)

# --------------------------------------------------------
# multiple_species_synthetic.snk

# snakemake -s multiple_species_synthetic.snk --configfile configs/test_multiple_synthetic.json --use-conda --conda-prefix {conda_prefix} --cores {cores} -d res/test_multiple_synthetic --rerun-incomplete -p --rerun-triggers mtime --scheduler greedy --notemp --keep-going -n


wd = os.path.dirname(os.getcwd())

jinfo = {
	'MGYG000000002': 10,
	'MGYG000000099': 20,
	'MGYG000000146': 24,
	'MGYG000005233': 17,
	'MGYG000015138': 15,
	'MGYG000025996': 6,
	'MGYG000037464': 7,
	'MGYG000133141': 14,
	'MGYG000258696': 14,
}

jdata = {
	'samples': {
		"uhgg": {
			sample: {
				'circular': 'true',
				'nanopore_qual': 'average',
				'quantity': coverage,
				'refpath': os.path.join(wd, f'dataset/uhgg/strains/{sample}.fa.gz'),
				'seed': 50
			}
			for sample, coverage in jinfo.items()
		}
	},

	'outputs': [
		{
			'group': 'all',
			'assembly': {
				'name': 'kraken_ref',
				'mode': 'nanopore'
			},
			'vcalling': 'longshot',
			'phasing': {
				'name': 'floria',
				'readtype': 'nanopore',
				'fmode': 'none',
				'post_assembler': 'wtdbg2',
				'assembler_rtype': 'long_reads',
				'assembler_preset': 'nano'
			},
			'ref_filtering':{
				'length': 'none',
				'similarity': 'none'
			}
		},

		{
			'group': 'all',
			'assembly': {
				'name': 'kraken_ref',
				'mode': 'nanopore'
			},
			'vcalling': 'longshot',
			'phasing': {
				'name': 'strainberry',
				'readtype': 'nanopore'
			},
			'ref_filtering':{
				'length': 'none',
				'similarity': 'none'
			}
		},


		{
			"group": "all",
			"assembly": {
				"name": "strainxpress",
				"mode": "fast"
			},
			"ref_filtering": {
				"length": "none",
				"similarity": "none"
			}
		},

		{
		"group": "all",
			"assembly": {
				"name": "megahit"
			},
			"ref_filtering": {
				"length": "none",
				"similarity": "none"
			}
		},

		{
		"group": "all",
			"assembly": {
				"name": "flye",
				"read": "nanopore",
				"mode": "raw"
			},
			"ref_filtering": {
				"length": "1000",
				"similarity": "99.5"
			},
			"refcomp": [
				'fastani', 'mummer'
			]
		}
	]
}

outfile = './test_multiple_synthetic.json'
with open(outfile, 'w') as f:
    json.dump(jdata, f, indent=4, sort_keys=True)