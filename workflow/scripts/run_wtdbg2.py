# @Author: jsgounot
# @Date:   2022-12-08 09:42:12
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-12-13 16:26:47

# Run both wtdbg2 and wtpoa
# Touch empty fasta file if input file or intermediate file are empty
# Compress the working directory to avoid having too much output file

import sys, os, shutil
from pathlib import Path

fastq   = snakemake.input.fastq
wdir    = snakemake.params.wdir
subname = snakemake.params.subname
preset  = snakemake.params.preset
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

if not isempty(fastq):
	cmdline = f'wtdbg2 {preset} -t {threads} -i {fastq} -o {of} 2> {logfile}'	
	os.system(cmdline)

	# usual empty gz file does not have st_size == 0 but it looks like
	# wtdbg2 just touch empty file with gz extension sometimes

	if not isempty(temp1):
		cmdline = f'wtpoa-cns -t {threads} -i {temp1} -fo {temp2} 2>> {logfile}'
		os.system(cmdline)
	else:
		Path(temp2).touch()

else:
	Path(temp2).touch()

outfile = snakemake.output.fa
archive = snakemake.output.ar

shutil.copyfile(temp2, outfile)

cmdline = f'tar -czf {archive} {wdir} --remove-files'
os.system(cmdline)

if rmfq == 'True':
	os.remove(fastq)