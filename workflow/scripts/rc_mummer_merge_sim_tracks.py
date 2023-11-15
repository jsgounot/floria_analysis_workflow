# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-02-17 16:48:10
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-08-11 14:07:43

'''
The goal of this script is to extract region with high similarity between the reference
genomes, defined as subtracks. It will also define the opposite tracks of these regions, defined
as revtracks.
'''

import os
import pandas as pd
from pybedtools import BedTool, Interval
from Bio import SeqIO

fname = snakemake.input[0]
sim = snakemake.wildcards['sim']
sim = float(sim[:2] + '.' + sim[2:])
length = int(snakemake.wildcards['length'])

simtracks = snakemake.output['simtracks']
revtracks = snakemake.output['revtracks']

def clean_filter(sdf):
    sdf.columns = ['contig', 'start', 'end']
    sdf['start'], sdf['end'] = sdf[['start', 'end']].min(axis=1), sdf[['start', 'end']].max(axis=1)
    return sdf
    
def bedtoolformdf(df, minlen=0, minidy=0):
    sdf = df[df['%IDY'] >= minidy]
    rdf = sdf[sdf['R_HITLEN'] >= minlen][['R_NAME', 'R_BEG', 'R_END']]
    rdf = clean_filter(rdf)
    qdf = sdf[sdf['Q_HITLEN'] >= minlen][['Q_NAME', 'Q_BEG', 'Q_END']]
    qdf = clean_filter(qdf)
    sdf = pd.concat((rdf, qdf)).sort_values(['contig', 'start', 'end'])
    return BedTool.from_dataframe(sdf)

print (f'Start extraction with sim: {sim} and length: {length}')

df = pd.read_csv(fname, sep='\t')
btdf = bedtoolformdf(df, length, sim)
btdf = btdf.merge()

ldf = len(btdf)
tc = btdf.total_coverage()
print (f'Extract {ldf:,} similarity tracks accounting for a total size of {tc:,} bp')
btdf.to_dataframe().to_csv(simtracks, sep='\t', index=False)

# Extracting opposite regions
refs = snakemake.params['refs']
refs = pd.DataFrame([
    {'refpath': os.path.basename(refpath), 'chrom': record.id, 
        'start': 0, 'end': len(record.seq)}
    for refpath in refs
    for record in SeqIO.parse(refpath, 'fasta')
    ])

refs = refs[['chrom', 'start', 'end', 'refpath']]
ldf = BedTool.from_dataframe(refs)
sdf = ldf.subtract(btdf)

ldf = len(sdf)
tc = sdf.total_coverage()
print (f'Reverse extract {ldf:,} non-similar tracks accounting for a total size of {tc:,} bp')
sdf.to_dataframe().to_csv(revtracks, sep='\t', index=False)