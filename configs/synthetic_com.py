# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2024-01-18 15:30:28
# @Last Modified by:   jsgounot
# @Last Modified time: 2024-03-11 17:00:16

import os
import json

scoverage = {
	'MGYG000000002': 10, 'MGYG000000038': 24, 'MGYG000000077': 16, 'MGYG000000099': 20, 
	'MGYG000000142': 22, 'MGYG000000146': 24, 'MGYG000000195': 14, 'MGYG000000196': 18, 
	'MGYG000000212': 23, 'MGYG000000262': 13, 'MGYG000001292': 17, 'MGYG000001300': 11, 
	'MGYG000001315': 17, 'MGYG000001338': 15, 'MGYG000001359': 21, 'MGYG000001378': 11, 
	'MGYG000002272': 14, 'MGYG000002279': 6, 'MGYG000002395': 5, 'MGYG000002396': 24, 
	'MGYG000002438': 10, 'MGYG000002469': 23, 'MGYG000002478': 22, 'MGYG000002492': 15, 
	'MGYG000002506': 24, 'MGYG000002545': 10, 'MGYG000003452': 8, 'MGYG000003683': 12, 
	'MGYG000003899': 10, 'MGYG000005127': 12, 'MGYG000005165': 7, 'MGYG000005233': 17, 
	'MGYG000005316': 5, 'MGYG000005570': 8, 'MGYG000005594': 5, 'MGYG000006035': 7, 
	'MGYG000006106': 10, 'MGYG000006146': 9, 'MGYG000006426': 16, 'MGYG000006533': 5, 
	'MGYG000006554': 16, 'MGYG000006657': 13, 'MGYG000006741': 5, 'MGYG000006778': 16, 
	'MGYG000006853': 7, 'MGYG000006888': 11, 'MGYG000006932': 10, 'MGYG000007036': 9, 
	'MGYG000010975': 18, 'MGYG000011485': 23, 'MGYG000011921': 22, 'MGYG000013885': 24, 
	'MGYG000015138': 15, 'MGYG000015151': 10, 'MGYG000018553': 8, 'MGYG000020051': 24, 
	'MGYG000021014': 12, 'MGYG000025754': 23, 'MGYG000025996': 6, 'MGYG000029665': 23, 
	'MGYG000032808': 6, 'MGYG000035951': 14, 'MGYG000037464': 7, 'MGYG000038970': 17, 
	'MGYG000041261': 6, 'MGYG000042722': 10, 'MGYG000050633': 21, 'MGYG000053833': 8, 
	'MGYG000059649': 14, 'MGYG000060819': 23, 'MGYG000067694': 10, 'MGYG000067967': 8, 
	'MGYG000068727': 13, 'MGYG000070171': 23, 'MGYG000071756': 12, 'MGYG000072213': 20, 
	'MGYG000077694': 11, 'MGYG000078292': 22, 'MGYG000082710': 21, 'MGYG000086726': 15, 
	'MGYG000087716': 6, 'MGYG000103793': 13, 'MGYG000112172': 9, 'MGYG000112749': 11, 
	'MGYG000120645': 9, 'MGYG000127504': 20, 'MGYG000132626': 13, 'MGYG000132680': 16, 
	'MGYG000133141': 14, 'MGYG000136675': 22, 'MGYG000142592': 13, 'MGYG000143066': 8, 
	'MGYG000146642': 24, 'MGYG000146993': 11, 'MGYG000152080': 9, 'MGYG000156166': 20, 
	'MGYG000161032': 5, 'MGYG000164526': 7, 'MGYG000164573': 15, 'MGYG000166196': 14, 
	'MGYG000170250': 14, 'MGYG000179292': 18, 'MGYG000185651': 24, 'MGYG000187110': 17, 
	'MGYG000195261': 19, 'MGYG000204845': 23, 'MGYG000204911': 19, 'MGYG000213752': 7, 
	'MGYG000216666': 5, 'MGYG000225022': 13, 'MGYG000225426': 17, 'MGYG000232699': 23, 
	'MGYG000245554': 5, 'MGYG000252297': 8, 'MGYG000257988': 20, 'MGYG000258696': 14, 
	'MGYG000260357': 10, 'MGYG000269216': 11, 'MGYG000273626': 22, 'MGYG000277300': 13
	}

wd = os.path.dirname(os.getcwd())

jdata = {
	'samples': {
		"uhgg1": {
			sample: {
				'circular': 'true',
				'nanopore_qual': 'average',
				'quantity': coverage,
				'refpath': os.path.join(wd, f'dataset/uhgg/strains/{sample}.fa.gz'),
				'seed': 50
			}
			for sample, coverage in scoverage.items()
		}
	},

	'outputs': [
		
		# Floria and Strainberry phasing
		# with Kraken as reference inputs

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
				'name': 'kraken_ref',
				'mode': 'illumina'
			},
			'vcalling': 'longshot',
			'phasing': {
				'name': 'floria',
				'readtype': 'illumina',
				'fmode': 'none',
				'post_assembler': 'megahit',
				'assembler_rtype': 'short_reads',
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

		# nanopore Floria and Strainberry phasing
		# with naive assemblies as reference

		{
			'group': 'all',
			'assembly': {
				'name': 'flye',
				'read': 'nanopore',
				'mode': 'raw'
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
				'name': 'megahit'
			},
			'vcalling': 'longshot',
			'phasing': {
				'name': 'floria',
				'readtype': 'illumina',
				'fmode': 'none',
				'post_assembler': 'megahit',
				'assembler_rtype': 'short_reads',
				'assembler_preset': 'none'
			},
			'ref_filtering':{
				'length': 'none',
				'similarity': 'none'
			}
		},
		

		# StrainXPress and naive assemblies

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
				"length": "none",
				"similarity": "none"
			},
			#"refcomp": [
			#	'fastani', 'mummer'
			#]
		}
	],

	'miscs': {
		'GENERIC_THREADS': 48
	}
}

outfile = './synthetic_com.json'
with open(outfile, 'w') as f:
    json.dump(jdata, f, indent=4, sort_keys=True)

