# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-04-25 11:47:00
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-06-16 11:44:39

# Need to check if I should something special for illumina reads
# for example if I need to conserve / remove other paired-read during cleaning

import pysam
import pandas as pd

fname = snakemake.input['maxrank']
df = pd.read_csv(fname, sep='\t', compression='gzip', low_memory=False)

def rename_contig(contig):
    contig = str(contig)
    if contig == 'nan': return 'Unmapped'
    return contig.split('|')[0].split('_')[0]

df['co'] = df['contig'].apply(rename_contig)
rc = df.groupby('co')['name'].apply(set).to_dict()
empty = set()

# --------------------------------------------------------

used = total = 0

fname = snakemake.input['bam']
outfile = snakemake.output[0]

with pysam.AlignmentFile(fname, "rb") as f:
    with pysam.AlignmentFile(outfile, 'wb', template=f) as o:
        for aln in f.fetch(until_eof=True):
            name = aln.query_name
            contig = rename_contig(aln.reference_name)
            if name in rc.get(contig, empty):
                o.write(aln)
                used += 1
            total += 1

prc = used * 100 / total

log = snakemake.log[0]
with open(log, 'w') as f:
    f.write(f'Used: {used} - Total: {total} - Prc: {prc:.2f}%')