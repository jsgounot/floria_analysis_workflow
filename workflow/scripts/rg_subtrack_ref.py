# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-02-24 18:35:37
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-03-06 18:04:44

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

fname = snakemake.input['reffasta']
#fname = 'res/synthetic_com_2_2/references/used/uhgg1/sample/MGYG000032808.fasta'
fdata = {record.id: record.seq for record in SeqIO.parse(fname, 'fasta')}
sname = os.path.basename(fname)

fname = snakemake.input['simtrack']
#fname = 'res/synthetic_com_2_2/refcomp/uhgg1.simtracks.reverse.995.1000.tsv.gz'
df = pd.read_csv(fname, sep='\t')
df = df[df['name'] == sname].sort_values(['chrom', 'start', 'end'])
df = df.reset_index(drop=True)

if len(df) == 0:
    raise Exception('Not implemented yet: One reference is completly overlapped by similar segments')

def make_record(idx, row, fdata):
    contig, start, end = row['chrom'], row['start'], row['end']
    seq = fdata[contig][start:end]
    rid = contig + '_' + str(idx)
    description = f'subtract {start} {end}'
    return SeqRecord(id=rid, seq=seq, description=description)

subdata = (make_record(idx, row, fdata) for idx, row in df.iterrows())
outfile = snakemake.output[0]
with open(outfile, 'w') as f:
    SeqIO.write(subdata, f, 'fasta')