# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-03-25 15:49:46
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-03-25 16:47:11

incfile = outfile = snakemake.input[0]
sample = snakemake.params[0] + '-'
outfile = snakemake.output[0]

import gzip

with gzip.open(incfile, 'rt') as fi:
	with gzip.open(outfile, 'wt') as fo:
		for idx, line in enumerate(fi):
			if idx % 4 == 0: 
				assert line.startswith('@')
				line = '@' + sample + line[1:]
			fo.write(line)