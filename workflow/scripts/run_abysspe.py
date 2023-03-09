# @Author: jsgounot
# @Date:   2022-12-08 09:42:12
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-01-09 13:00:01

import sys, os, shutil
from pathlib import Path

r1      = snakemake.input.r1
r2      = snakemake.input.r2
wdir    = snakemake.params.wdir
rmfq    = snakemake.params.rmfq
logfile = snakemake.log
threads = snakemake.threads

isempty = lambda fname: os.path.isfile(fname) == False or os.stat(fname).st_size == 0

name = os.path.join(wdir, 'assembly')
scaffold = name + '-scaffolds.fa'

if os.path.isdir(wdir):
	shutil.rmtree(wdir) 

os.makedirs(wdir)

if not (isempty(r1) and isempty(r2)):
	cmdline = f'abyss-pe name={name} k=96 B=2G in=\'{r1} {r2}\' --quiet 2> {logfile}'
	os.system(cmdline)

	# sometime you have a segmentation dump 
	# I guess when there are too few reads
	try: scaffold = os.readlink(scaffold)
	except FileNotFoundError: Path(scaffold).touch()	

else:
	Path(scaffold).touch()

outfile = snakemake.output.fa
archive = snakemake.output.ar

try :
    shutil.copyfile(scaffold, outfile)
except BlockingIOError:
    # Weird error I have on the cluster
    cmdline = f'cp {scaffold} {outfile}'
    os.system(cmdline)

cmdline = f'tar -czf {archive} {wdir} --remove-files'
os.system(cmdline)

if rmfq == 'True':
	os.remove(r1)
	os.remove(r2)
