# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-04-04 10:11:10
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-04-04 15:07:37

import glob, os
from pathlib import Path

dname = os.path.dirname
bname = os.path.basename

fnames = 'phasing/glopp/*/*/*/haplostats.tsv'
fnames = glob.glob(fnames)

for fname in fnames:
	group = bname(dname(dname(dname(fname)))) 
	mode  = bname(dname(dname(fname)))
	read  = bname(dname(fname))

	fname = Path(fname).absolute()

	outfile = f'reports/glopp/haplotigs/{group}.{mode}.{read}.tsv'
	if not os.path.isfile(outfile): os.symlink(fname, outfile)

# ----------------------------------------------------------------------

fnames = 'stats/assemblies/*/*/mummer/report.tsv'
fnames = glob.glob(fnames)

for fname in fnames:
	group = bname(dname(dname(dname(fname)))) 
	fid = bname(dname(dname(fname)))
	outfile = f'reports/assemblies/{group}.{fid}.tsv'

	fname = Path(fname).absolute()

	if not os.path.isfile(outfile): os.symlink(fname, outfile)

# ----------------------------------------------------------------------

fnames = 'stats/assemblies/*/*/mummer/circos'
fnames = glob.glob(fnames)

for fname in fnames:
	group = bname(dname(dname(dname(fname)))) 
	fid = bname(dname(dname(fname)))
	outfile = f'reports/circos/{group}.{fid}'

	fname = Path(fname).absolute()
	if not os.path.exists(outfile): os.symlink(fname, outfile)