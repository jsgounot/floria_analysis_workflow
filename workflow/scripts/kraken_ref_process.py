# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-04-24 10:28:01
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-06-16 11:14:19

# note regarding kraken readlength formatting between illumina and nanopore
# using the --paired option (my guess), readlength in the kraken output file (not the report)
# has this form '{readlength read1}|{readlength read2}'.

import gzip
import pandas as pd

fname = snakemake.input['meta']
meta = pd.read_csv(fname, sep='\t', index_col=0)
tlen = meta.groupby('taxid')['seqlen'].sum()

fname = snakemake.input['report']
names = ['abundance', '#covered', '#assigned', 'rank', 'taxid', 'name']
report = pd.read_csv(fname, sep='\t', names=names)
report['name'] = report['name'].str.strip()

fname = snakemake.input['output']
names = ['classified', 'readid', 'taxid', 'size']
df = pd.read_csv(fname, sep='\t', names=names, compression='gzip', usecols=[0,1,2,3])

# illumina formatting read length (see script file head)
if snakemake.wildcards.krtype == 'illumina':
	df['size'] = df['size'].apply(lambda value: sum(int(s) for s in value.split('|')) )

df = df.groupby('taxid')['size'].agg(('size', 'sum'))
df.columns = ['assigned', 'readlength']
df = df.reset_index()
df['assigned_abu'] =  (df['assigned'] * 100 / df['assigned'].sum()).round(2)

df['reflen'] = df['taxid'].map(tlen)
df['cov'] = df['readlength'] / df['reflen']

report = report.merge(df, on='taxid', how='left')
report['reflen_mb'] =  report['reflen'] / 1e6

outfile = snakemake.output[0]
report.to_csv(outfile, sep='\t')

mincov = float(snakemake.params['mincov'])
tids = (str(tid) for tid in df[df['cov'] >= mincov]['taxid'])
outfile = snakemake.output[1]
with open(outfile, 'w') as f:
	f.write('\n'.join(tids))
