# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-06-19 16:46:40
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-06-19 16:49:18

import gzip
import pandas as pd
from Bio import SeqIO

data = []
fname = snakemake.input[0]
with gzip.open(fname, 'rt') as f:
	fdata = SeqIO.parse(f, 'fasta')
	for record in fdata:
		seqname, _, taxid = record.id.split('|')
		seqlen = len(record.seq)
		data.append({
			'seqname': seqname, 'taxid': taxid, 'seqlen': seqlen
			})

df = pd.DataFrame(data)
outfile = snakemake.output[0]
df.to_csv(outfile, sep='\t', compression='gzip')