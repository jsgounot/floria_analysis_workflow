# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-04-18 15:55:04
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-04-18 15:57:10

from Bio import SeqIO

fname = snakemake.input[0]
outfile = snakemake.output[0]
sample = snakemake.params['sample']
sample = sample.replace(' ', '_')

fdata = list(SeqIO.parse(fname, 'fasta'))
for record in fdata:
	record.id = record.id + '_' + sample

with open(outfile, 'w') as f:
	SeqIO.write(fdata, f, 'fasta')