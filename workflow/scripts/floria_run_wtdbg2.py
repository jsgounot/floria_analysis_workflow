# @Author: jsgounot
# @Date:   2022-12-08 09:42:12
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-06-26 18:24:25

# Run both wtdbg2 and wtpoa
# Touch empty fasta file if input file or intermediate file are empty
# Compress the working directory to avoid having too much output file

import sys, os, shutil, glob, gzip, tarfile, time
import concurrent.futures, subprocess

import logging
from logging.handlers import RotatingFileHandler

from pathlib import Path
from Bio import SeqIO

bname = os.path.basename
dname = os.path.dirname

# Logger
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s')
 
logpath = snakemake.log[0]
file_handler = RotatingFileHandler(logpath, 'w')
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
 
#stream_handler = logging.StreamHandler()
#stream_handler.setLevel(logging.DEBUG)
#logger.addHandler(stream_handler)

isempty = lambda fname: os.path.isfile(fname) == False or os.stat(fname).st_size == 0

def make_haplotype_assembly(elements):
    tar, fastq, preset, wthreads = elements

    fastq = fastq.replace('.fastq.gz', '.wtdbg2.fastq.gz')
    assert os.path.isfile(fastq)
    logger.info(f'Run wtdbg2 with: {fastq}')

    wdir = fastq[:-9]
    subname = os.path.basename(fastq)[:-9]
    sub_outfile = wdir + '.fasta'
    archive = wdir + '.tar.gz'
    outfile_cp = sub_outfile + '.gz'
    logfile = os.path.join(wdir, 'wtdbg2.logfile.txt')

    if os.path.isfile(outfile_cp):
        logger.warning(f'File {outfile_cp} already processed (rerun?), continue ...')
        return outfile_cp

    if not os.path.isfile(fastq):
        logger.error(f'Fastq file not extracted? {fastq}')
        raise Exception('Fastq file not extracted? {fastq}')

    of = os.path.join(wdir, subname)
    temp1 = of + '.ctg.lay.gz'
    temp2 = of + '.ctg.fa'

    if os.path.isdir(wdir):
        shutil.rmtree(wdir) 

    os.makedirs(wdir)

    with open(logfile, 'w') as f:

        if not isempty(fastq):
            cmdline = f'wtdbg2 {preset} -t {wthreads} -i "{fastq}" -o "{of}"'
            subprocess.call(cmdline, shell=True, stdout=f, stderr=f)

            # usual empty gz file does not have st_size == 0 but it looks like
            # wtdbg2 just touch empty file with gz extension sometimes

            if not isempty(temp1):
                cmdline = f'wtpoa-cns -t {wthreads} -i "{temp1}" -fo "{temp2}"'
                subprocess.call(cmdline, shell=True, stdout=f, stderr=f)

            else:
                Path(temp2).touch()

        else:
            Path(temp2).touch()

    shutil.copyfile(temp2, sub_outfile)
    shutil.rmtree(wdir)
    os.remove(fastq)

    cmdline = f'gzip -f "{sub_outfile}"'
    os.system(cmdline)

    return outfile_cp

def iter_contigs(fnames):
    for fname in fnames:
        logger.info(f'Read contig file: {fname}')
        contig = bname(dname(dname(fname)))
        part = bname(fname).split('.')[0][:-5] # we remove the '_part'
        
        with gzip.open(fname, 'rt') as f:
            fdata = SeqIO.parse(f, 'fasta')
            for record in fdata:
                record.id = contig + '__' + part + '__' + record.id
                yield record

        os.remove(fname)

# ------------------------------------------------------------------------------------------
# Tar content parsing

reads_tar = snakemake.input['reads_ar']
outfile = snakemake.output[0]
threads = snakemake.threads
preset  = snakemake.params['preset']

subprocess_thread = min(threads, 16)
workers = threads // subprocess_thread
workers = workers or 1

logger.info(f'Run wtdbg2 using reads from: {reads_tar}')
logger.info(f'Final contig file path: {outfile}')
logger.info(f'Threads: {threads} - Workers: {workers}')
logger.info(f'WTDBG2 preset options: {preset}')

rtype = snakemake.wildcards['rtype']
validate = lambda fname: fname.endswith('.fastq.gz') and bname(dname(fname)) == rtype

tar = tarfile.open(reads_tar)
fnames = [(reads_tar, path, preset, subprocess_thread) for path in tar.getnames() if validate(path)]

assert fnames

# ------------------------------------------------------------------------------------------
# Concurrent processing of all parts

results = []

with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
    for resfname in executor.map(make_haplotype_assembly, fnames):
        results.append(resfname)

if os.path.isfile(outfile):
    os.remove(outfile)

records = iter_contigs(results)
with gzip.open(outfile, 'wt') as f:
    SeqIO.write(records, f, 'fasta')
