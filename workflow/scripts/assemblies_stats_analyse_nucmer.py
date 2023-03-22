# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-03-28 19:32:34
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-12-08 14:26:47


'''
About this:
We try here to assign one query contig to its best reference
This is based on this strainberry script:
https://github.com/rvicedomini/strainberry-analyses/blob/fcc370d61293f05b2c09b18bb2130c1f11985ad2/workflow/scripts/assembly_stats.py

For each query / reference, a score is calculated based on the average identity x the number of based aligned
At least 50% of the reference sequence / genome must be covered

Original contains an issue. Here I use bedtool and pandas, the code is a bit slow but we don't really care in this case.
'''

import os, glob

bname = os.path.basename
dname = os.path.dirname

from pybedtools import BedTool
import pandas as pd

def load_qcoords(fname):
    names = ['R_BEG', 'R_END', 'Q_BEG', 'Q_END', 'R_HITLEN', 'Q_HITLEN', 
            '%IDY', 'R_LEN', 'Q_LEN', 'R_COV', 'Q_COV', 'R_NAME', 'Q_NAME']
    
    df = pd.read_csv(fname, sep='\t', names=names)
    df['fidx'] = bname(fname).split('.')[0]
    return df

def get_coverage(df):
    if len(df) > 1:
        # quite slow but does the job correctly
        # https://github.com/rvicedomini/strainberry-analyses/issues/1
        bt = BedTool.from_dataframe(df[['Q_NAME', 'Q_MIN', 'Q_MAX']].sort_values('Q_MIN'))
        return bt.total_coverage()
    else:
        return (df['Q_MAX'] - df['Q_MIN']).sum()

def extract_sim(df):
    df['Q_MIN'] = df[['Q_BEG', 'Q_END']].min(axis=1) - 1
    df['Q_MAX'] = df[['Q_BEG', 'Q_END']].max(axis=1)

    df['weight'] = df['R_HITLEN'] + df['Q_HITLEN']

    keys = ['Q_NAME', 'fidx']
    cov = df.groupby(keys).apply(get_coverage)
    cov = cov.rename('aligned_bases').reset_index()

    df['widy'] = df['%IDY'] * df['weight']
    avgIdy = df.groupby(keys)['widy'].sum() / df.groupby(keys)['weight'].sum()
    avgIdy = avgIdy.rename('avg_idy').reset_index()

    qlen = df.groupby('Q_NAME')['Q_LEN'].nunique()
    assert qlen[qlen != 1].empty
    qlen = df.drop_duplicates('Q_NAME').set_index('Q_NAME')['Q_LEN'].to_dict()

    sdf = avgIdy.merge(cov, how='outer')
    sdf['size_query'] = sdf['Q_NAME'].map(qlen)

    sdf['prc_query'] = (100 * sdf['aligned_bases']) / sdf['size_query']

    return sdf

def filter_best(df):
    sdf = extract_sim(df)

    # Minimum 50% of the reference bases covered
    sdf = sdf[sdf['prc_query'] > 50]

    # We select the best alignment based on a score
    # score = average identity * aligned bases
    # are we sure of this?
    sdf['score'] = sdf['avg_idy'] * sdf['aligned_bases']
    sdf['best'] = sdf.groupby('Q_NAME')['score'].transform(max) == sdf['score']

    # What should we do with duplicates?
    # I keep them here but in original strainberry script, \
    # they are removed randomly and only one is kept
    # sdf = sdf[sdf['best']].drop_duplicates('Q_NAME')
    return sdf

def main(fnames, outfile):
    fnames = sorted(fnames)
    df = pd.concat([load_qcoords(fname) for fname in fnames])
    df = filter_best(df)
    df.to_csv(outfile, sep='\t', index=False)  

fnames = snakemake.input
outfile = snakemake.output[0]
main(fnames, outfile)


