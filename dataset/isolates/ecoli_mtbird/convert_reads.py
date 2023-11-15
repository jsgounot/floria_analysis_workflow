# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-04-21 09:06:38
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-11-10 10:08:21

import os
import gzip
import glob

import concurrent.futures

'''
Reads are not named correctly for bwa or minimap2 to recognize them as paired
In what I downloaded: @SRR11361693.1.1 1 length=101
The fixed version: @SRR11361693.1 1 length=101/1
'''

def rename_reads(fname, outfile):
    strain = os.path.basename(os.path.dirname(fname))
    
    ridx = fname[-10]
    assert ridx in '12'
    ridx = '/' + ridx + '\n'
    
    with gzip.open(fname, 'rt') as f:
        with gzip.open(outfile, 'wt') as of:
            while True:

                try: line = next(f)
                except StopIteration : break

                line = line.strip().split()
                line = line[0][:-1] + strain + ridx
                of.write(line)

                line = next(f)
                of.write(line)

                line = next(f)
                line = line.strip().split()
                line = line[0][:-1] + strain + ridx
                of.write(line)

                line = next(f)
                of.write(line)

    print (f'Done with {fname}')

fname = snakemake.input[0]
outfile = snakemake.output[0]
rename_reads(fname, outfile)