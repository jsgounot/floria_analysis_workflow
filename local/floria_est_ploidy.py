# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-06-22 16:06:15
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-06-22 16:24:43

import glob, os
import pandas as pd

df = []

directory = 'synthetic_com'
fnames = os.path.join('res', directory, 'phasing/floria/presplit/*/*/*/wdir/contig_ploidy_info.tsv')
fnames = glob.glob(fnames)

for fname in fnames:
	split = fname.split('/')
	sdf = pd.read_csv(fname, sep='\t')
	sdf['cat'] = 'presplit'
	sdf['group'] = split[5]
	sdf['seqtype'] = split[6]
	sdf['vcaller'] = split[7]
	sdf['tid'] = 'None'
	df.append(sdf)

fnames = os.path.join('res', directory, 'phasing/floria/split/*/*/*/*/wdir/contig_ploidy_info.tsv')
fnames = glob.glob(fnames)

for fname in fnames:
	split = fname.split('/')
	sdf = pd.read_csv(fname, sep='\t')
	sdf['cat'] = 'split'
	sdf['group'] = split[5]
	sdf['seqtype'] = split[6]
	sdf['vcaller'] = split[8]
	sdf['tid'] = split[7]
	df.append(sdf)

df = pd.concat(df)
outfile = os.path.join('res', directory, 'reports', 'floria_est_ploidy.tsv')
df.to_csv(outfile, sep='\t')