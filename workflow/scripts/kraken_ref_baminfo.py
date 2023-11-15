# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-04-24 13:04:54
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-07-24 14:48:39

import pysam
import pandas as pd
import numpy as np

# pd.set_option('display.max_colwidth', None)

CIGAR = 'MIDSNHPCXBE'

def iter_bam_contig(bamfile):
    # Massively drop memory consumption to chunk by contigs

    bamfile = pysam.AlignmentFile(bamfile, "rb")
    rows = []
    tracked_contig = None

    for aln in bamfile.fetch(until_eof=True):
        iql = aln.infer_query_length() or np.nan
        irl = aln.infer_read_length() or np.nan
        len_ops, num_ops = aln.get_cigar_stats()

        mapq = aln.mapping_quality,
        flag = aln.flag,
        assert len(mapq) == 1
        assert len(flag) == 1
        mapq = mapq[0]
        flag = flag[0]

        cigar = {
            CIGAR[idx]: value
            for idx, value in enumerate(len_ops)
        }

        name = aln.query_name
        contig = aln.reference_name
        taxid = int(contig.split('_')[-1]) if contig else 'None'

        if tracked_contig != contig and rows:
            #print (f'Process {tracked_contig} {len(rows)}')
            yield pd.DataFrame(rows)
            rows = []
            tracked_contig = contig

        rows.append({
            'name': name,
            'contig': contig,
            'taxid': taxid,
            'iql': iql,
            'irl': irl,
            'mapq': mapq,
            'flag': flag,
            ** cigar
            })

    if rows:
        yield pd.DataFrame(rows)

    bamfile.close()

# First pass
bamfile = snakemake.input[0]

df = pd.concat(iter_bam_contig(bamfile))
df = df.reset_index(drop=True)
df['sim'] = (df['M'] / df['irl']).fillna(0) * 100

outfile = snakemake.output['baminfo']
df.to_csv(outfile, sep='\t', compression='gzip', index=False)

# -------- 

df = df[df['sim'] == df.groupby(['name', 'taxid'])['sim'].transform(max)]
df['rank'] = df.groupby('name')['sim'].rank(ascending=False, method='dense')

sdf = df[df['rank'] == 1][['name', 'contig', 'taxid']]
outfile = snakemake.output['maxrank']
sdf.to_csv(outfile, sep='\t', compression='gzip', index=False)

# -------- 

df = df.sort_values('name')

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