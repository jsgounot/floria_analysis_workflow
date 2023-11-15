# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-07-17 12:37:56
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-07-17 14:29:46

import pandas as pd
import gzip, os

pd.set_option('display.max_colwidth', None)

bname = os.path.basename
dname = os.path.dirname

def read_fastq(fdata):
    readid = None
    for idx, line in enumerate(fdata):
        if readid:
            line = line.decode("utf-8").strip()
            yield {'rid': readid, 'rlen': len(line)}
            readid = None
        if idx % 4 == 0:
            readid = line.decode("utf-8").split()[0].strip()[1:]
            readid = readid[:-2] if readid[-2] == '/' else readid

def iter_reads_info_gzips(fnames, rids):
    for fname in fnames:
        bname = os.path.basename(fname)
        with gzip.open(fname) as f:
            for readinfo in read_fastq(f):
                if readinfo['rid'] in rids: continue
                yield readinfo
                # yield {'bname': bname,  ** readinfo}

def get_annot_cinfos(fname):
    sdf = pd.read_csv(fname, sep='\t', compression='gzip')
    sdf['tid'] = bname(dname(dname(fname)))
    return sdf

cinfos = snakemake.input['cinfos']
cinfos = pd.concat((get_annot_cinfos(fname) for fname in cinfos))
rids = set(cinfos['rid'])

reads = snakemake.input['reads']
allreads = pd.DataFrame(iter_reads_info_gzips(reads, rids))
allreads =  allreads.drop_duplicates('rid')
allreads['haplo_id'] = 'missing_reads'

allreads = pd.concat((cinfos, allreads))
outfile = snakemake.output[0]
allreads.to_csv(outfile, sep='\t', compression='gzip', index=False)