# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-03-28 15:58:48
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-11-30 10:52:52

'''
I first wanted to include this at the end of assemblies_stats_analyse_nucmer.py script,
but this does not work well with snakemake to manage unknown number of output.
'''

import pandas as pd
from Bio import SeqIO

def main(table, query, idx, outfile):
    df = pd.read_csv(table, sep='\t')

    sdf = df[(df['fidx'].astype(str) == str(idx)) & (df['best'])]
    ids = set(sdf['Q_NAME'].astype(str))

    # In some cases we can have no ids for an idx
    # This means we did not find a single contig
    # which align the best to this reference
    # in this case we just write an empty file
    if not ids:
        with open(outfile, 'w'):
            pass
        return

    fdata = SeqIO.parse(query, 'fasta')
    with open(outfile, 'w') as f:
        frec = sorted((record for record in fdata if record.id in ids),
            key = lambda rec: rec.id)

        if not len(frec) == len(ids):
            raise Exception(f'Error: Unable to fetch all contig found in nucmer results ({len(frec)} vs {len(ids)})')

        SeqIO.write(frec, f, 'fasta')

table = snakemake.input[0]
query = snakemake.input[1]
idx = snakemake.params[0]
outfile = snakemake.output[0]
main(table, query, idx, outfile)


