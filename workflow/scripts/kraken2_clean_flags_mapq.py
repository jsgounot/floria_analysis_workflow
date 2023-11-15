# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-07-18 15:10:19
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-07-18 15:27:57

used = total = 0
fname = snakemake.input['bam']
outfile = snakemake.output[0]

best_mapq = defaultdict(lambda: defaultdict(int))
with pysam.AlignmentFile(fname, "rb") as f:
    for aln in f.fetch(until_eof=True):
        name = aln.query_name
        contig = rename_contig(aln.reference_name)
        mapq = aln.mapping_quality
        value = best_mapq[contig][name]
        if mapq > value: best_mapq[contig][name] = mapq

best_mapq = {contig: {name: (mapq, True) for name, mapq in cdata.items()}
    for contig, cdata in best_mapq.items()}

with pysam.AlignmentFile(fname, "rb") as f:
    with pysam.AlignmentFile(outfile, 'wb', template=f) as o:
        for aln in f.fetch(until_eof=True):
            name = aln.query_name
            contig = rename_contig(aln.reference_name)
            mapq = aln.mapping_quality
            
            if name in rc.get(contig, empty):

                # we check from primary alignment
                rmapq, singleton = best_mapq[contig][name]
                if mapq == rmapq and singleton:
                    aln.is_secondary = False
                    best_mapq[contig][name] = (mapq, False)
                else:
                    aln.is_secondary = True

                o.write(aln)
                used += 1
            total += 1

prc = used * 100 / total

log = snakemake.log[0]
with open(log, 'w') as f:
    f.write(f'Used: {used} - Total: {total} - Prc: {prc:.2f}%')