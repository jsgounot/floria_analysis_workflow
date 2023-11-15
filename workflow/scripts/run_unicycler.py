# @Author: jsgounot
# @Date:   2022-12-08 09:42:12
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-01-09 12:39:52

import sys, os, shutil
from pathlib import Path

r1      = snakemake.input.r1
r2      = snakemake.input.r2
wdir    = snakemake.params.wdir
subname = snakemake.params.subname
rmfq    = snakemake.params.rmfq
logfile = snakemake.log
threads = snakemake.threads

isempty = lambda fname: os.path.isfile(fname) == False or os.stat(fname).st_size == 0

of = os.path.join(wdir, subname)
temp1 = of + '.ctg.lay.gz'
temp2 = of + '.ctg.fa'

if os.path.isdir(wdir):
	shutil.rmtree(wdir) 

os.makedirs(wdir)

if not (isempty(r1) and isempty(r2)):
	cmdline = f'unicycler --short1 {r1} --short2 {r2} --threads {threads} --out {wdir}'	
	os.system(cmdline)

else:
	Path(temp2).touch()

outfile = snakemake.output.fa
archive = snakemake.output.ar

shutil.copyfile(temp2, outfile)

cmdline = f'tar -czf {archive} {wdir} --remove-files'
os.system(cmdline)

if rmfq == 'True':
	os.remove(r1)
	os.remove(r2)