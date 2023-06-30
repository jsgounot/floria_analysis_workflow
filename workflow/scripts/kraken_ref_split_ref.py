# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-05-23 13:36:40
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-05-23 15:34:13

import gzip
from Bio import SeqIO

lc = snakemake.input['lc']
fa = snakemake.input['fa']
of = snakemake.output[0]

with open(lc) as f:
	contigs = set(line.strip() 
		for line in f)

with gzip.open(fa, 'rt') as f:
	fdata = SeqIO.parse(f, 'fasta')
	fdata = [record for record in fdata if record.id in contigs]
	assert len(fdata) == len(contigs)

SeqIO.write(fdata, of, 'fasta')