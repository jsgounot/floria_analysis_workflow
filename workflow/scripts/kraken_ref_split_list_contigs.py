# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-05-23 13:36:40
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-06-06 14:16:05

import os
from collections import defaultdict

dic = defaultdict(set)

fname = snakemake.input['lc']

with open(fname) as f:
	for line in f:
		line = line.strip()
		print (line)
		taxid = 'tid' + line.split('_')[-1]
		dic[taxid].add(line)

dname = snakemake.params['dname']
if not os.path.isdir(dname): os.mkdir(dname)

krtype = snakemake.wildcards.krtype

for taxid, contigs in dic.items():
	outfile = os.path.join(dname, f'{taxid}.{krtype}.contigs.txt')
	with open(outfile, 'w') as f:
		f.write('\n'.join(sorted(contigs)))

with open(snakemake.output[0], 'w') as f:
	f.write('\n'.join(sorted(set(dic))))