# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-07-18 16:54:53
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-07-18 17:09:58

from Bio import SeqIO
from collections import defaultdict
import gzip, os

fasta = snakemake.input[0]
dname = snakemake.params['dname']

dname = snakemake.params['dname']
if not os.path.isdir(dname): os.mkdir(dname)
krtype = snakemake.wildcards.krtype

fdata = defaultdict(list)
with gzip.open(fasta, 'rt') as f:
	for record in SeqIO.parse(f, 'fasta'):
		taxid = record.id.split('_')[-1]
		fdata[taxid].append(record)

for taxid, records in fdata.items():
	outfile = os.path.join(dname, f'{taxid}.{krtype}.fa.gz')
	with gzip.open(outfile, 'wt') as f:
		SeqIO.write(records, f, 'fasta')

with open(snakemake.output[0], 'w') as f:
	f.write('\n'.join(sorted(fdata)))