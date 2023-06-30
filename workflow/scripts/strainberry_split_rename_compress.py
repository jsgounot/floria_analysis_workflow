# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-02-17 16:48:10
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-02-17 16:53:26

from Bio import SeqIO
import gzip

fname = snakemake.input['scaffolds']
tid = snakemake.wildcards['tid']
outfile = snakemake.output[0]

fdata = list(SeqIO.parse(fname, 'fasta'))

for record in fdata:
	record.id = f'{tid}_{record.id}'

with gzip.open(outfile, 'wt') as f:
	SeqIO.write(fdata, f, 'fasta')