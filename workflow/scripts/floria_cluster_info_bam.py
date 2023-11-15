# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-07-17 10:20:37
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-07-17 13:52:06

import pysam
import pandas as pd
import tarfile, fnmatch, gzip, os

# pd.set_option('display.max_colwidth', None)

bname = os.path.basename
dname = os.path.dirname

fname = snakemake.input['archive']
rows = []

with tarfile.open(fname, "r:*") as tar:
    for fname in tar.getnames():
        if fname.endswith('.haplosets'):
            contig = bname(dname(fname))

            f = tar.extractfile(fname)
            for line in f:
                line = line.decode("utf-8").strip()
                if line.startswith('>'):
                    haplo_id = line.split('.')[0][1:]
                    continue
                else:
                    rid, start, end = line.split()
                    rows.append({
                        'contig': contig, 'haplo_id': haplo_id, 'rid': rid
                        })

df = pd.DataFrame(rows)

fname = snakemake.input['bam']
bamfile = pysam.AlignmentFile(fname, "rb")
allreads = {read.query_name: (read.reference_name, read.query_length or read.infer_query_length()) 
    for read in bamfile.fetch()}

# Note: Here contig assigned to the read is the last one observed in the bamfile!

missing = set(allreads) - set(df['rid'])
sdf = pd.DataFrame((
    {
        'contig': 'missing', 'haplo_id': 'missing_floria', 
        'rid': rid, 'contig': allreads[rid][0]
    }
    for rid in missing)
)

df = pd.concat((df, sdf))
df['rlen'] = df['rid'].apply(lambda rid: allreads[rid][1])

outfile = snakemake.output[0]
df.to_csv(outfile, sep='\t', index=False, compression='gzip')