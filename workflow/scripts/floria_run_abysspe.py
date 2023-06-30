# @Author: jsgounot
# @Date:   2022-12-08 09:42:12
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-06-26 18:25:45

import sys, os, shutil, glob, gzip, tarfile
import concurrent.futures, subprocess

import logging
from logging.handlers import RotatingFileHandler

from collections import defaultdict
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

def make_haplotype_assembly(elements, wthreads=1):

    tar, partname = elements

    logger.info(f'Work with: {partname}')

    wdir = partname + '.abysspe'
    subname = os.path.basename(partname)
    outfile = wdir + '.fasta'
    archive = wdir + '.tar.gz'
    outfile_cp = outfile + '.gz'

    if os.path.isfile(outfile_cp):
        f'File {outfile_cp} already processed (rerun?), continue ...'
        return outfile_cp

    r1 = partname + '_paired1.abysspe.fastq.gz'
    r2 = partname + '_paired2.abysspe.fastq.gz'

    '''
    # Similar to wtdbg2, the decompression is done with snakemake now 

    tar = tarfile.open(tar)
    for member in tar.getmembers():
        if member.name == r1 or member.name == r2:
            tar.extract(member)

    '''

    assert os.path.isfile(r1)
    assert os.path.isfile(r2)

    name = os.path.join(wdir, 'assembly')
    scaffold = name + '-scaffolds.fa'
    logfile = os.path.join(wdir, 'abysspe.logfile.txt')

    if os.path.isdir(wdir):
        shutil.rmtree(wdir) 

    os.makedirs(wdir)

    # Update 05.05.2023 - After discussion with Jim
    # I change k=96 to k=48. The initial value was too important for some
    # non simulated reads, leading to almost no assembly at all with glopp.
    # see: https://github.com/bcgsc/abyss#assembling-using-a-paired-de-bruijn-graph

    if not (isempty(r1) and isempty(r2)):
        cmdline = f'abyss-pe name={name} k=48 B=2G in=\'{r1} {r2}\' j={wthreads} --quiet'

        with open(logfile, 'w') as f:
            subprocess.call(cmdline, shell=True, stdout=f, stderr=f)

        # sometime you have a segmentation dump 
        # I guess when there are too few reads
        try: scaffold = os.readlink(scaffold)
        except FileNotFoundError: Path(scaffold).touch()    

    else:
        Path(scaffold).touch()

    try :
        shutil.copyfile(scaffold, outfile)
    except BlockingIOError:
        # Weird error I have on clusters
        cmdline = f'cp {scaffold} {outfile}'
        os.system(cmdline)

    #cmdline = f'tar -czf {archive} {wdir} --remove-files'
    #os.system(cmdline)

    shutil.rmtree(wdir) 
    os.remove(r1)
    os.remove(r2)

    cmdline = f'gzip -f "{outfile}"'
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

# ------------------------------------------------------------------------------------------
# Tar content parsing

reads_tar = snakemake.input['reads_ar']
outfile = snakemake.output[0]
threads = snakemake.threads

logger.info(f'Run abyss-pe using reads from: {reads_tar}')
logger.info(f'Final contig file path: {outfile}')
logger.info(f'Threads: {threads}')

tar = tarfile.open(reads_tar)
tar_files = tar.getmembers()
parts = defaultdict(dict)

rtype = snakemake.wildcards['rtype']
validate = lambda fname: fname.endswith('.fastq.gz') and bname(dname(fname)) == rtype

for path in tar.getnames():
    if validate(path):
        part = path[:-17]
        rid = path[-10]

        if rid == '1': parts[part]['r1'] = path
        elif rid == '2': parts[part]['r2'] = path
        else: raise Exception('Unable to find rid from path?')

assert parts

for part, value in parts.items():
    assert 'r1' in value and 'r2' in value

nparts = len(parts)
logger.info(f'Process {nparts} parts')

parts = ((reads_tar, part) for part, value in parts.items())

# ------------------------------------------------------------------------------------------
# Concurrent processing of all parts

results = []
workers = threads

with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
    for resfname in executor.map(make_haplotype_assembly, parts):
        results.append(resfname)

if os.path.isfile(outfile):
    os.remove(outfile)

records = iter_contigs(results)
with gzip.open(outfile, 'wt') as f:
    SeqIO.write(records, f, 'fasta')
