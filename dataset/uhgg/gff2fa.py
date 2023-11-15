# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-04-24 10:46:51
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-04-27 17:46:22

import gzip

fname = snakemake.input[0]
outfi = snakemake.output[0]

with gzip.open(fname, 'rt') as f:
	record = False
	with gzip.open(outfi, 'wt') as fo:
		for line in f:
			if record:
				fo.write(line)
			if line.startswith('##FASTA'):
				record = True

assert record == True