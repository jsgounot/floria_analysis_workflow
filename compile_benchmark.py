# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-01-17 17:05:14
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-01-17 17:32:38

import glob
import pandas as pd
import os

flist = [
	'res/*/benchmarks/*',
	'res/*/benchmarks/*/*',
	'res/*/benchmarks/*/*/*'
]

values = []

for idx, fnames in enumerate(flist):
	fnames = glob.glob(fnames)

	for fname in fnames:
		print (fname)

		if not os.path.isfile(fname):
			continue

		bname = os.path.basename(fname)
		group = fname.split('/')[1]

		sdf = pd.read_csv(fname, sep='\t')
		sdf['group'] = group
		sdf['bname'] = bname
		sdf['subgroup'] = os.path.basename(os.path.dirname(fname)) if idx == 2 else 'main'
		values.append(sdf)

df = pd.concat(values)
for group, sdf in df.groupby('group'):
	outfile = f'res/{group}/reports/benchmark.tsv'
	sdf.to_csv(outfile, sep='\t')