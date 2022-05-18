# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-05-18 15:29:35
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-05-18 15:53:00

import glob
import os

fnames = snakemake.params[0]
suffix = snakemake.params[1]
preset = snakemake.params[2]
core_j = snakemake.params[3]
outfile = snakemake.output[0]


with open(outfile, 'w') as f:
	fnames = glob.glob(fnames)
	for fname in fnames:

		if not os.path.isfile(fname):
			raise Exception(f'File does not exist: {fname}')

		if os.stat(fname).st_size==0:
			# Sometime we have empty files, we juste ignore them
			continue

		outdir = fname[:-6] + '_' + suffix
		os.makedirs(outdir, exist_ok=True)

		bname = os.path.basename(fname[:-6]) + '.' + suffix
		outfile = os.path.join(outdir, bname)

		cmdline = f'(wtdbg2 {preset} -t {core_j} -i {fname} -o {outfile} && wtpoa-cns -t {core_j} -i {outfile}.ctg.lay.gz -fo {outfile}.ctg.fa) 2> {outfile}.log.txt \n'
		f.write(cmdline)
