# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-03-28 10:44:06
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-03-28 16:00:42

import os, glob

bname = os.path.basename
dname = os.path.dirname

from Bio import SeqIO

def iter_contigs(fnames):
	fnames = glob.glob(fnames)
	for fname in fnames:
		print ('Read contig file:', fname)
		contig = bname(dname(dname(dname(fname))))
		part = bname(fname).split('_')[0]
		fdata = SeqIO.parse(fname, 'fasta')
		for record in fdata:
			record.id = contig + '__' + part + '__' + record.id
			yield record

incfile = snakemake.input[0]
fnames = snakemake.params[0]
outfile = snakemake.output[0]

records = iter_contigs(fnames)
with open(outfile, 'w') as f:
	SeqIO.write(records, f, 'fasta')