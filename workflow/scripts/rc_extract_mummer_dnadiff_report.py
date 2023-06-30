# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-02-17 16:48:10
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-02-17 16:53:26

import os
import pandas as pd

fnames = snakemake.input['dnadiff']
outfile = snakemake.output[0]
refs = snakemake.params['refs']

def getprc(report_value):
    left_index = report_value.index('(') + 1
    return float(report_value[left_index:-2])

mummerres = []

for fname in fnames:
    with open(fname) as f:
        idx1, idx2, _ = os.path.basename(fname).split('.')
        subres = {'idx1': idx1, 'idx2': idx2,}
        for line in f:
            if line.startswith('AlignedBases'):
                line = line.strip().split()
                prc1, prc2 = getprc(line[1]), getprc(line[2])
                subres['prc_aligned1'], subres['prc_aligned2'] = prc1, prc2
                
            elif line.startswith('AvgIdentity'):
                line = line.strip().split()
                prc1, prc2 = float(line[1]), float(line[2])
                subres['avg_identity1'], subres['avg_identity2'] = prc1, prc2
        
        mummerres.append(subres)
                
ab = pd.DataFrame(mummerres)
ab.sort_values(['idx1', 'idx2']).head()

refs = {str(idx): fname for idx, fname in enumerate(refs)}
ab['ref1'] = ab['idx1'].map(refs)
ab['ref2'] = ab['idx2'].map(refs)

ab.to_csv(outfile, sep='\t')