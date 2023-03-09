# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-03-28 10:44:06
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-12-08 11:34:31

import os, glob

bname = os.path.basename
dname = os.path.dirname

from Bio import SeqIO

def iter_contigs(fnames):
	for fname in fnames:
		print ('Read contig file:', fname)
		contig = bname(dname(dname(fname)))
		part = bname(fname).split('.')[1][:-5] # we remove the '_part'
		fdata = SeqIO.parse(fname, 'fasta')
		for record in fdata:
			record.id = contig + '__' + part + '__' + record.id
			yield record

fnames = snakemake.input
outfile = snakemake.output[0]

records = iter_contigs(fnames)
with open(outfile, 'w') as f:
	SeqIO.write(records, f, 'fasta')