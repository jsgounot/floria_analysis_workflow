# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-03-28 15:43:51
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-04-04 11:04:57

import itertools
import gzip
import random
import shutil

def read_file(fname):
    if fname.endswith('.gz'):
        with gzip.open(fname, 'rt') as f:
            yield from f
    else:
        with open(fname) as f:
            yield from f

def extract_read_sizes(fname):
    rsizes = []
    record = False
    for idx, line in enumerate(read_file(fname)):
        if record:
            rsizes.append(len(line.strip()))
            record = False

        if idx % 4 == 0:
            line = line.strip()
            record = True

    return rsizes

def write_outfile(readfile, outfile, use_idxs):
    with gzip.open(outfile, 'wb') as f:
        for idx, line in enumerate(read_file(readfile)):
            if idx % 4 == 0:
                if idx // 4 in use_idxs:
                    record = True
                else:
                    record = False

            if record:
                f.write(line.encode())

ref   = snakemake.input['ref']
r1    = snakemake.input['r1']
r2    = snakemake.input['r2']
cov   = snakemake.params['cov']
seed  = snakemake.params['seed']

if cov == 'all':
    r1o = snakemake.output['r1o']
    shutil.copyfile(reads[0], r1o)

    r2o = snakemake.output['r1o']
    shutil.copyfile(reads[1], r2o)

    exit(0)

cov = int(cov)
random.seed(seed)

refsize = sum(len(line.strip()) for line in read_file(ref) if not line.startswith('>'))
target_len = refsize * cov

rsizes_one = extract_read_sizes(r1)
rsizes_two = extract_read_sizes(r2)
rsizes = [first + rsizes_two[idx] for idx, first in enumerate(rsizes_one)]

rrange = [idx for idx in range(0, len(rsizes))]
random.shuffle(rrange)

cumsize = 0
for ridx, idx in enumerate(rrange):
    cumsize += rsizes[idx]
    if cumsize >= target_len:
        break

use_idxs = set(rrange[:ridx + 1])
tcov = cumsize / refsize

print (f'Use of {len(use_idxs)} paired-reads out of {len(rsizes)}')
print (f'Mean coverage: {tcov:.1f}  - ref size: {refsize} - reads cumulative length: {cumsize:,}')

if cumsize < target_len:
    msg = f'Not enough reads to hit target coverage ({cov}X, {target_len:,}, {ref})'
    raise Exception(msg)

r1o = snakemake.output['r1o']
write_outfile(r1, r1o, use_idxs)

r2o = snakemake.output['r2o']
write_outfile(r2, r2o, use_idxs)
