# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-04-24 13:04:54
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-06-21 16:49:11

import pysam
import pandas as pd
import numpy as np

# First pass
bamfile = snakemake.input[0]
bamfile = pysam.AlignmentFile(bamfile, "rb")

CIGAR = 'MIDSNHPCXBE'
df = []

for read in bamfile.fetch(until_eof=True):
    iql = read.infer_query_length()
    irl = read.infer_read_length()
    len_ops, num_ops = read.get_cigar_stats()

    cigar = {
    	CIGAR[idx]: value
    	for idx, value in enumerate(len_ops)
    }

    name = read.query_name
    contig = read.reference_name
    taxid = int(contig.split('_')[-1]) if contig else 'None'

    df.append({
    	'name': name,
    	'contig': contig,
    	'taxid': taxid,
    	'iql': iql,
    	'irl': irl,
    	** cigar
    	})

bamfile.close()

df = pd.DataFrame(df)
outfile = snakemake.output['baminfo']
df.to_csv(outfile, sep='\t', compression='gzip', index=False)

# -------- 

df['sim'] = (df['M'] / df['irl']).fillna(0) * 100
df = df[df['sim'] == df.groupby(['name', 'taxid'])['sim'].transform(max)]

# fun = lambda serie: serie.rank(ascending=False)
# df['rank'] = df.groupby('name')['sim'].transform(fun)

df['rank'] = df.groupby('name')['sim'].rank(ascending=False, method='first')

outfile = snakemake.output['maxrank']
df.to_csv(outfile, sep='\t', compression='gzip', index=False)

# -------- 

df = df.sort_values('name')

'''
Old code, should produce the same results but much slower. I keep it if
something looks wrong with the new way.

def transform_name(sdf):
    sdf['rank'] = sdf['sim'].rank(ascending=False, method='first')
    sdf = sdf[(sdf['rank'] == 1) | (sdf['rank'] == 2)].sort_values('rank')
    sim2 = sdf.iloc[1]['sim'] if len(sdf) > 1 else np.nan
    serie = sdf.iloc[0]
    serie['sim_rank2'] = sim2
    return serie

df = df.groupby('name').apply(transform_name)
df = df.reset_index(drop=True)

outfile = snakemake.output['filtered']
df.to_csv(outfile, sep='\t', compression='gzip', index=False)
'''

sdf = df[df['rank'] < 3]
sdf = sdf.drop_duplicates(['name', 'rank'])
sdf['contig'] = sdf['contig'].fillna('Unmapped')

sim_df = sdf.pivot_table(index='name', columns='rank', values='sim', fill_value=0).reset_index()

if len(sim_df.columns == 2):
    # Happen if we don't have any co-alignment
    sim_df[2] = 0

sim_df.columns = ['name', 'sim_1', 'sim_2']

fun = lambda values: values
co_df = sdf.pivot_table(index='name', columns='rank', values='contig', aggfunc=fun).reset_index()

if len(co_df.columns == 2):
    # Happen if we don't have any co-alignment
    co_df[2] = np.nan

co_df.columns = ['name', 'contig_1', 'contig_2']

sdf = sim_df.merge(co_df, on='name')

outfile = snakemake.output['top2']
sdf.to_csv(outfile, sep='\t', compression='gzip', index=False)