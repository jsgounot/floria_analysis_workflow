# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-02-17 16:48:10
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-02-17 16:53:26

import os
import pandas as pd

fnames = snakemake.input['coors']
outfile = snakemake.output[0]
refs = snakemake.params['refs']

names = ['R_BEG', 'R_END', 'Q_BEG', 'Q_END', 'R_HITLEN', 'Q_HITLEN', 
            '%IDY', 'R_LEN', 'Q_LEN', 'FRM', 'TAGS', 'R_NAME', 'Q_NAME']

mummerres = []

def iter_df(fnames):
    for fname in fnames:
        with open(fname) as f:
            idx1, idx2, _ = os.path.basename(fname).split('.')
            df = pd.read_csv(fname, sep='\t', names=names, skiprows=4)
            df['idx1'], df['idx2'] = idx1, idx2
            mummerres.append(df)
            yield df
                
ab = pd.concat(iter_df(fnames))
ab.sort_values(['idx1', 'idx2']).head()

refs = {str(idx): fname for idx, fname in enumerate(refs)}
ab['ref1'] = ab['idx1'].map(refs)
ab['ref2'] = ab['idx2'].map(refs)

ab.to_csv(outfile, sep='\t', index=False)