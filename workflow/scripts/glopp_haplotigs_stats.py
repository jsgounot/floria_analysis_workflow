# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-03-26 09:55:04
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-04-04 10:06:23

import glob, os, sys
import pandas as pd
from collections import Counter

SFUNS = [
    lambda rname: '_'.join(rname.split('_')[:2]),
    lambda rname: rname.split('-')[0]
    ]

def read_all_parts(fname):
    hid = None
    d = {}

    with open(fname) as f:
        for line in f:
            if line.startswith('#'):
                line = line[1:].strip().split(',')
                hid, coverage, error_rate = line
                d[hid] = {
                    'coverage': coverage,
                    'error_rate': error_rate,
                    'reads': []
                }

            else:
                assert hid is not None
                line = line.strip().split('\t')
                d[hid]['reads'].append(line)

    return d

def reads_extract_sample(hdata, sfun):
    for read in hdata['reads']:
        #print (read[0], sfun(read[0]))
        yield sfun(read[0])

def reads_positions(hdata):
    if len(hdata['reads']) == 0:
        yield -1

    for read in hdata['reads']:
        for position in read[1:]:
            yield int(position)
            
def make_parts_stats(ap, sfun):
    d = []
    samples = set()

    for haplotig, hdata in ap.items():
        d.append({
            'hid': haplotig,
            'coverage': hdata['coverage'],
            'error_rate': hdata['error_rate'],
            'min_pos': min(reads_positions(hdata)),
            'max_pos': max(reads_positions(hdata)),
            }
            )

        counter = Counter(reads_extract_sample(hdata, sfun))
        samples |= set(counter)

        assert len(counter) < 10
        d[-1].update(counter)

    df = pd.DataFrame(d)

    sdf = df[list(samples)]
    df['max_prop'] = sdf.max(axis=1) / sdf.sum(axis=1)
    df['min_prop'] = sdf.min(axis=1) / sdf.sum(axis=1)

    return df

def extract_haplofile_info(fnames):
    d = []
    for fname in fnames:
        print (fname)
        with open(fname) as f:
            line = next(f).strip()
            hid, supp = line.split(',')
            start, stop = supp.split('.')

            positions = [int(line.strip().split('\t')[0].split(':')[1])
            for line in f]

            if len(positions) == 0: positions.append(-1)

            fpos, lpos = min(positions), max(positions)

            d.append({
                'hid': hid[1:], 'hap_first_snp': start, 'hap_last_snp': stop,
                'hap_first_pos': fpos, 'hap_last_pos': lpos
                })

    return pd.DataFrame(d)

def main(dname, outfile, sfun):
    fnames = glob.glob(os.path.join(dname, '*', 'all_part.txt'))
    results = []

    for fname in fnames:
        contig = os.path.basename(os.path.dirname(fname))
        ap = read_all_parts(fname)

        stat = make_parts_stats(ap, sfun)
        stat['contig'] = contig

        subfnames = os.path.join(os.path.dirname(fname), 'haplotypes', '*_hap.txt')
        subfnames = glob.glob(subfnames)
        subfnames = extract_haplofile_info(subfnames)

        stat = stat.merge(subfnames, on='hid', how='outer')
        results.append(stat)

    df = pd.concat(results)

    print (df)

    corder = ['contig', 'hid', 'coverage', 'error_rate', 'min_pos', 'max_pos', 'hap_first_snp', 'hap_last_snp', 'hap_first_pos', 'hap_last_pos', 'min_prop', 'max_prop']
    corder = corder + sorted(set(df.columns) - set(corder))

    df = df[corder]
    df.to_csv(outfile, sep='\t')

    print (df)

dname = sys.argv[1]
spatt = sys.argv[2]
sfun = SFUNS[int(spatt)]

outfi = os.path.join(dname, 'haplostats.tsv')
main(dname, outfi, sfun)