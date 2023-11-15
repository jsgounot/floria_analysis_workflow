# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-04-12 17:39:07
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-04-12 17:44:44

from Bio import SeqIO

fname = snakemake.input[0]
outfile = snakemake.output[0]

fdata = SeqIO.parse(fname, 'fasta')
maxs = sums = 0
mainc = None

for record in fdata:
	lseq = len(record.seq)
	sums += lseq
	if lseq > maxs:
		maxs = lseq
		mainc = record

lseq = len(mainc.seq)
prc   = lseq * 100 / sums
print (f'{lseq} - {prc:.2f}')
assert prc > 90

fdata = [mainc]
with open(outfile, 'w') as f:
	SeqIO.write(fdata, f, 'fasta')