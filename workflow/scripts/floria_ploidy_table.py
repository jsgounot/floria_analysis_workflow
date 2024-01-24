# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-06-28 13:38:56
# @Last Modified by:   jsgounot
# @Last Modified time: 2024-01-22 13:19:33

import tarfile, os
import pandas as pd
from collections import defaultdict
from Bio import SeqIO

tar_arc = snakemake.input['archive']
fasta = snakemake.input['fasta']
kraken = snakemake.input['kraken']
depth = snakemake.input['depth']
vcf  = snakemake.input['vcf']

df = None
with tarfile.open(tar_arc, "r:*") as tar: 
	for fname in tar.getnames():
		if os.path.basename(fname) == 'contig_ploidy_info.tsv':
			if df is not None: raise Exception('Multiple contig ploidy files')
			df = pd.read_csv(tar.extractfile(fname), sep='\t')

assert df is not None

# Contig size
flen = {record.id: len(record.seq) for record in SeqIO.parse(fasta, 'fasta')}
df['contig_size'] = df['contig'].map(flen)

# Kraken results
df['taxid'] = df['contig'].apply(lambda name: int(name.split('_')[-1]))
kdf = pd.read_csv(kraken, sep='\t', index_col=0)
kdf.columns = ['kraken_' + column if column != 'taxid' else column for column in kdf.columns]
df = df.merge(kdf, on='taxid', how='left')

# Depth
ddf = pd.read_csv(depth, sep='\t', compression='gzip', index_col=0)
ddf = ddf.groupby('contig')[['npos', 'sum_depth']].sum().reset_index()
ddf['mean_cov'] = ddf['sum_depth'] / ddf['npos']
ddf.columns = ['bdepth_' + column if column != 'contig' else column for column in ddf.columns]
df = df.merge(ddf, on='contig', how='left')

# VCF
d_confident = defaultdict(list)
d_ambiguous = defaultdict(list)

with open(vcf, 'r') as f:
	for line in f:
		if line.startswith('#'): continue
		line = line.strip().split('\t')
		contig, info = line[0], line[7]
		info = dict(e.split('=') for e in info.split(';') if e)
		info['AC'] = info['AC'].split(',')[1]

		DP, AC, AM = int(info['DP']), int(info['AC']), int(info['AM'])
		d_confident[contig].append(AC / (DP - AM))
		d_ambiguous[contig].append(AC / DP)

d_nbase = {contig: len(cvalues) for contig, cvalues in d_confident.items()}
d_confident = {contig: sum(cvalues) / len(cvalues) for contig, cvalues in d_confident.items()}
d_ambiguous = {contig: sum(cvalues) / len(cvalues) for contig, cvalues in d_ambiguous.items()}

df['af_nba'] = df['contig'].map(d_nbase)
df['af_con'] = df['contig'].map(d_confident)
df['af_amb'] = df['contig'].map(d_ambiguous)

tid = snakemake.params.get('tid', None)
if tid: df['tid'] = tid

outfile = snakemake.output[0]
df.to_csv(outfile, sep='\t', compression='gzip', index=False)