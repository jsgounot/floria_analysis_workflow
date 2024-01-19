# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2024-01-17 14:19:21
# @Last Modified by:   jsgounot
# @Last Modified time: 2024-01-18 13:55:12

# ----------------------------------------------------------------------
# Synthetic pairwise divergence

# 7 strains based on Strainberry strains used for Figure 4 (variable divergence panel)
# are used here. They are going to be downloaded and reads simulated automaticaly with the pipeline
# for both Floria and Strainberry

import json
from itertools import combinations

strains = {
	'K_12': 'NC_000913.3',
	'ME8067': 'NZ_CP028703.1',
	'LD27_1': 'NZ_CP047594.1', 
	'Y5': 'NZ_CP013483.1', 
	'EC590': 'NZ_CP016182.2', 
	'RM14721': 'NZ_CP027105.1', 
	'AMSCJX03': 'NZ_CP058355.1', 
	'H5': 'NZ_CP010169.1'
}

jdata = {
	'samples': {
		strain: {
			'K_12': {
				'ncbinuc': strains['K_12'],
				'mapping_ref': 'True',
				'quantity': 20,
				'circular': 'true',
				'nanopore_qual': 'average',
				'seed': 50
			},

			strain: {
				'ncbinuc': identifier,
				'mapping_ref': 'False',
				'quantity': 30,
				'circular': 'true',
				'nanopore_qual': 'average',
				'seed': 50
			}
		}

		for strain, identifier in strains.items()
		if strain != 'K_12'

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
				'post_assembler': 'flye',
				'assembler_rtype': 'long_reads',
				'assembler_preset': 'none'
			},
			'ref_filtering':{
				'length': 'none',
				'similarity': 'none'
			},
			'refcomp': [
				'fastani', 'mummer'
				]
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
		},

		{
			'group': 'all',
			'assembly': {
				'name': 'flye',
				'read': 'nanopore',
				'mode': 'raw'
			},
			'ref_filtering':{
				'length': 'none',
				'similarity': 'none'
			}
		}
	],

	'miscs': {
		'GENERIC_THREADS': 4
	}
}

outfile = './ecoli_divergence.json'
with open(outfile, 'w') as f:
    json.dump(jdata, f, indent=4, sort_keys=True)