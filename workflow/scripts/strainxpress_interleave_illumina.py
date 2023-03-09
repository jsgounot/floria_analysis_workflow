# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-04-05 15:41:06
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-05-17 15:01:53

# Copy from here: https://github.com/ekg/interleave-fastq/blob/master/interleave-fastq

import sys, os

def interleave(f1, f2, of):
    # Interleaves two (open) fastq files.
    while True:
        line = f1.readline()
        if line.strip() == "":
            break

        of.write(line.strip() + '\n')
        
        for i in range(3):
            of.write(f1.readline().strip() + '\n')
        
        for i in range(4):
            of.write(f2.readline().strip() + '\n')


file1 = snakemake.input[0]
file2 = snakemake.input[1]
outfile = snakemake.output[0]

if file1[-2:] == "gz":
    import gzip
    with gzip.open(file1, 'rt') as f1:
        with gzip.open(file2, 'rt') as f2:
            with open(outfile, 'w') as of:
                interleave(f1, f2, of)

else:
    with open(file1) as f1:
        with open(file2) as f2:
            with open(outfile, 'w') as of:
                interleave(f1, f2, of)

