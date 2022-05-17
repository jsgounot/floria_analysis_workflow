# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-04-05 15:41:06
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-05-09 17:16:06

# Copy from here: https://github.com/ekg/interleave-fastq/blob/master/interleave-fastq

import sys, os

def interleave(f1, f2):
    # Interleaves two (open) fastq files.
    while True:
        line = f1.readline()
        if line.strip() == "":
            break

        print (line.strip())
        
        for i in range(3):
            print (f1.readline().strip())
        
        for i in range(4):
            print (f2.readline().strip())

if __name__ == '__main__':
    try:
        file1 = sys.argv[1]
        file2 = sys.argv[2]
    except:
        print ('Unable to parse arguments (python strainxpress_interleave_illumina.py fastq1 fastq2)')
        sys.exit(1)

    isfile = os.path.isfile
    if not isfile(file1) or not isfile(file2):
        raise OSError('Files not found')

    if file1[-2:] == "gz":
        import gzip
        with gzip.open(file1, 'rt') as f1:
            with gzip.open(file2, 'rt') as f2:
                interleave(f1, f2)

    else:
        with open(file1) as f1:
            with open(file2) as f2:
                interleave(f1, f2)

