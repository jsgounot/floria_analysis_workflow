# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-03-25 15:49:46
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-05-22 14:16:04

incfile = outfile = snakemake.input[0]
sample = snakemake.params[0] + '-'
outfile = snakemake.output[0]

import gzip

with gzip.open(incfile, 'rt') as fi:
	with gzip.open(outfile, 'wt') as fo:
		for idx, line in enumerate(fi):
			if idx % 4 == 0: 
				if not line.startswith('@'):
					raise Exception(f'Line does not start with @ at index {idx}')
				line = '@' + sample + line[1:]
			fo.write(line)