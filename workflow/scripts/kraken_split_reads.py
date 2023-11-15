# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-07-20 15:59:11
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-07-25 14:11:24

import gzip

readnames = snakemake.input['readnames']
reads = snakemake.input['reads']
outfile = snakemake.output[0]

with open(readnames) as f:
	names = {line.strip() for line in f}

found = set()
with gzip.open(reads, 'rt') as f:
	with gzip.open(outfile, 'wt') as of:
		record = False
		for idx, line in enumerate(f):
			if idx % 4 == 0:
				sline = line.strip().split()[0]
				sline = sline[1:-2] if sline[-2] == '/' else sline[1:]
				if sline in names:
					found.add(sline)
					record = True
				else:
					record = False

			if record:
				of.write(line)

if len(names) != len(found):
	raise Exception(f'Found {len(found)} records out of {len(names)}!')
