import os
from Bio import SeqIO

import pandas as pd

def parse_dnadiff_report(fname):
    res = {}
    with open(fname) as f:
        for line in f:
            line = line.strip().split()
            if line and line[0] not in res:
                res[line[0]] = line
    return res


def dnadiff_dup_bases(fname):
    count = 0
    with open(fname) as f:
        for line in f:
            line = line.strip().split('\t')
            if line[1] == 'DUP':
                count += int(line[4])
    return count

def fasta_n50(fname, gsize=0):
    fdata = SeqIO.parse(fname, 'fasta')
    lengths = sorted((len(record.seq) for record in fdata), reverse=True)
    mid = (gsize or sum(length)) / 2

    cumsum = 0
    for length in lengths:
        cumsum += length
        if cumsum > mid :
            return length

def parse_mummer_outputs(report, reflist):
    prefix = report[:-7]
    result = {}

    index = int(os.path.basename(report).split('.')[0])
    refpath = reflist[index]
    refname = os.path.splitext(os.path.basename(refpath))[0]

    report = parse_dnadiff_report(prefix + '.report')

    result['refidx'] = index
    result['refpath'] = refpath
    result['refbname'] = refname

    result['seq_num']  = int(report['TotalSeqs'][2])
    result['ref_size'] = int(report['TotalBases'][1])
    result['asm_size'] = int(report['TotalBases'][2])
    result['n50']      = fasta_n50(prefix + '.fa', result['ref_size'])

    result['unaligned_ref_bases'] = int(report['UnalignedBases'][1].split('(')[0])
    result['unaligned_ref'] = 100.0 * result['unaligned_ref_bases'] / result['ref_size']
    result['unaligned_asm_bases'] = int(report['UnalignedBases'][2].split('(')[0])
    result['unaligned_asm'] = 100.0 * result['unaligned_asm_bases'] / result['asm_size']

    result['ani'] = float(report['AvgIdentity'][1])
    
    result['aligned_ref_bases'] = int(report['AlignedBases'][1].split('(')[0])
    result['aligned_asm_bases'] = int(report['AlignedBases'][2].split('(')[0])
    result['dup_ratio'] = result['aligned_asm_bases'] / result['aligned_ref_bases']
    
    result['dup_bases'] = dnadiff_dup_bases(prefix + '.qdiff')
    result['cmp_bases'] = dnadiff_dup_bases(prefix + '.rdiff')
    
    result['snps']        = int(report['TotalSNPs'][1])
    result['inversions']  = int(report['Inversions'][2])
    result['relocations'] = int(report['Relocations'][2])
    result['transloc']    = int(report['Translocations'][2])
    
    return result

def main(reports, reflist, outfile):
    df = pd.DataFrame([
        parse_mummer_outputs(report, reflist)
        for report in reports
        ])

    df.to_csv(outfile, sep='\t', index=False)

reports = snakemake.input
reflist = snakemake.params[0]
outfile  = snakemake.output[0]

main(reports, reflist, outfile)