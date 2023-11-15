# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-07-24 10:37:42
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-07-24 11:05:13

# Replace any unusual base to N.
# Some variant callers do not like unusual bases (like freebayes)

import sys
from Bio import SeqIO
from Bio.Seq import Seq

def filter_clean(fname):
    fdata = SeqIO.parse(fname, 'fasta')

    usual = set('ATGCN')
    ident = set()

    for record in fdata:
        seq = str(record.seq).upper()
        found = set(seq)
        ident |= found

        for unusual in found - usual:
            seq = seq.replace(unusual, 'N')

        record.seq = Seq(seq)
        yield record
    
    # print (ident - usual)

fname = snakemake.input[0]
outfile = snakemake.output[0]

with open(outfile, 'w') as of: 
    SeqIO.write(filter_clean(fname), of, 'fasta')