# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-04-25 11:47:00
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-07-18 16:17:22

# Need to check if I should something special for illumina reads
# for example if I need to conserve / remove other paired-read during cleaning

import pysam
import pandas as pd

pd.set_option('display.max_colwidth', None)

fname = snakemake.input['maxrank']
df = pd.read_csv(fname, sep='\t', compression='gzip', low_memory=False)

def rename_contig(contig):
    contig = str(contig)
    if contig == 'nan': return 'Unmapped'
    return contig.split('|')[0].split('_')[0]


# ALL best values are kept. This means a read can be found in multiple taxid
# If you want to take one randomly, use method='first'
# If you want to ignore read which have multiple best sim, use average (rank will be >1 but <2)
df['rank'] = df.groupby('name')['sim'].rank(ascending=False, method='dense')
df = df[df['rank'] == 1]

df['co'] = df['contig'].apply(rename_contig)
rc = df.groupby('co')['name'].apply(set).to_dict()
empty = set()

# Important note
# This is not the cleanest way to do this. Since the mapq score is partially defined
# by other alignment, it is possible that the mapq score would be different (probably
# a bit better) for some alignments here if we remap reads against the subset of reference.
# However, most of the mapq will remain high as we assume here that most of the time the 
# best sim is also the primary alignment (by definition for minimap2 at least, the primary 
# aln is the one with the best mapq).
# print (df['mapq'].describe())

# It is very important to note that this bam file is not good enough when splitted
# since we use here a 'dense' version, we can have multiple reads mapping to different
# taxids. We need to clean bam file a second time to redefine primary and secondary.

test_name = 'MGYG000143066-218d32a7-d144-653a-f27e-01da784fb4f9'
test_name = 'MGYG000143066-3885a90c-94a7-28b9-be4e-878ca2caaeb1'

used = total = 0
fname = snakemake.input['bam']
outfile = snakemake.output[0]

best_mapq = {}
with pysam.AlignmentFile(fname, "rb") as f:
    for aln in f.fetch(until_eof=True):
        name = aln.query_name
        contig = rename_contig(aln.reference_name)
        mapq = aln.mapping_quality
        value = best_mapq.setdefault(name, 0)
        if mapq > value: best_mapq[name] = mapq

print (best_mapq[test_name])

with pysam.AlignmentFile(fname, "rb") as f:
    with pysam.AlignmentFile(outfile, 'wb', template=f) as o:
        for aln in f.fetch(until_eof=True):
            name = aln.query_name
            contig = rename_contig(aln.reference_name)
            if name in rc.get(contig, empty):
                mapq = aln.mapping_quality
                bmapq = best_mapq.get(name, None)

                if name == test_name:
                    print ('enter', aln.is_secondary, aln.flag)

                not_supp = aln.is_supplementary == False

                if not_supp and bmapq is not None and bmapq == mapq:
                    aln.is_secondary = False
                    del best_mapq[name]
                    if name == test_name:
                        print ('yes', aln.is_secondary, aln.flag)

                else:
                    aln.is_secondary = True
                    if name == test_name:
                        print ('no', aln.is_secondary, aln.flag)

                o.write(aln)
                used += 1
            total += 1

prc = used * 100 / total

log = snakemake.log[0]
with open(log, 'w') as f:
    f.write(f'Used: {used} - Total: {total} - Prc: {prc:.2f}%')
