# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-03-28 13:40:18
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-03-28 14:15:53

import pysam
import pysamstats
from scipy.stats import binom
import sys

bam_file_path = snakemake.input['aln']
ref_path =  snakemake.input['fasta']
outfile = snakemake.output[0]

binomial_null_error_rate = 0.04
p_value = 0.010
minimum_minor_al = 3

with pysam.AlignmentFile(bam_file_path) as bamfile:
    with open(outfile, 'w') as outf:

        outf.write('##fileformat=VCFv4.2\n')
        for contig in bamfile.references:
            outf.write(f'##contig=<{contig}>')
        
        tsize = sum(bamfile.lengths)

        outf.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth of reads passing MAPQ filter">\n')
        outf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n')
        for idx, record in enumerate(pysamstats.stat_variation(bamfile, fafile=ref_path)):

            if idx and idx % 100000 == 0:
                prc = idx * 100 / tsize
                print (f'{idx} {prc:.1f}%')

            contig = record['chrom']
            position = record['pos'] + 1
            ref = record['ref']

            altsc = {base: record[base] for base in 'ACGT' if base != ref}
            refc = record.get(ref, 0)

            probs = {alt: binom.cdf(count, refc + count, 0.02)
                for alt, count in altsc.items()
                if count >= minimum_minor_al}

            alts = [alt
                for alt, prob in probs.items()
                if 1 - prob < p_value
                ]

            if alts:
                dpcount = sum(altsc[alt] for alt in alts)
                outf.write(f'{contig}\t{position}\t.\t{ref}\t{",".join(alts)}\t1\tPASS\tDP={dpcount}\n')