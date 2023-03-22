# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-01-25 14:55:16
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-03-22 13:23:02

import glob, sys
import shutil

fnames = 'res/*/phasing/assemblies/*/*.glopp.*.fa.gz'
fnames = set(glob.glob(fnames))

dnames = 'res/*/phasing/glopp/*/inpref/*'
dnames = glob.glob(dnames)

try:
	rm = sys.argv[1] == 'rm'
except IndexError:
	rm = False

c = {True: 0, False: 0}

# inpref.glopp.hybrid_gp2763.wtdbg2.long_reads.longshot.fa.gz
# res/zb_D6331/phasing/glopp/zb_D6331_var/inpref/hybrid_gp765

for dname in dnames:
	dsplit = dname.split('/')
	maindir, subdir, gdir = dsplit[1], dsplit[4], dsplit[6]
	search = f'res/{maindir}/phasing/assemblies/{subdir}/inpref.glopp.{gdir}.wtdbg2.long_reads.longshot.fa.gz'
	c[search in fnames] += 1

	if search not in fnames:
		print (dname)

	if search not in fnames and rm:
		shutil.rmtree(dname)

print (c)