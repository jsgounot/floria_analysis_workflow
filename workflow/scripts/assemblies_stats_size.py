# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-02-28 17:53:46
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-02-28 17:57:21

import gzip, json
from Bio import SeqIO
from collections import defaultdict

fname = snakemake.input[0]
count = defaultdict(int)

with gzip.open(fname, 'rt') as f:
	records = SeqIO.parse(f, 'fasta')
	for record in records:
		count[len(record.seq)] += 1

outfile = snakemake.output[0]
with open(outfile, 'w') as f:
	json.dump(count, outfile, indent=2)