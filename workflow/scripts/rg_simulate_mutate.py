# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-02-24 18:35:37
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-03-06 18:04:44

from Bio import SeqIO
from Bio.Seq import Seq

import sys, random

random.seed(1)

def main(fasta, mutation_freq, outfile):

    records = list(SeqIO.parse(fasta, 'fasta'))

    for record in records:
        seq = list(record.seq)
        for idx, base in enumerate(record.seq):
            val = random.random()
            if val < mutation_freq:
                seq[idx] = random.choice([x for x in "ACTG" if x != base.upper()])
        record.seq = Seq(''.join(seq))

    with open(outfile, 'w') as f:
        SeqIO.write(records, f, 'fasta')

fname = snakemake.input[0]
mfreq = snakemake.params['snps_rate']

outfile = snakemake.output[0]
main(fname, mfreq, outfile)