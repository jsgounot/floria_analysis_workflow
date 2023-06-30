# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-01-19 15:35:44
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-01-19 15:44:03

import glob, os, sys
from itertools import combinations
from collections import defaultdict, Counter
from Bio import SeqIO

dname = sys.argv[1]
groups = os.path.join(dname, 'references/used/*/')
groups = sorted(glob.glob(groups))

for group in groups:
	print (f'{os.path.basename(os.path.dirname(group))}:')
	fnames = os.path.join(group, 'sample/*.fasta')
	fnames = glob.glob(fnames)

	# not really right but should be in most cases
	ref = fnames[0]

	sites = defaultdict(int)

	for a, b in combinations(fnames, 2):
		fa = {record.id: record.seq for record in SeqIO.parse(a, 'fasta')}
		fb = {record.id: record.seq for record in SeqIO.parse(b, 'fasta')}
		diff = tot = 0
		for recorid, aseq in fa.items():
			bseq = fb.get(recorid, None)
			if bseq and len(bseq) == len(aseq):
				for index, base in enumerate(aseq):
					if base != bseq[index]:
						diff += 1
						if a == ref: sites[index] += 1
					tot += 1

		if tot:
			prc = diff * 100 / tot
			print (f'{os.path.basename(a)} {os.path.basename(b)} {diff} {tot} {prc:.1f}')
		else:
			print (f'{os.path.basename(a)} {os.path.basename(b)} Nothing to compare')

	if sites:
		prc = len(sites) * 100 / tot
		print (f'# of polymorphic sites: {len(sites)}, prc: {prc:.1f}%')

		allelic = Counter(sites.values())
		print ('Allelic distribution:')
		for a, b in allelic.items():
			prc = b * 100 / sum(allelic.values())
			print (f'# Alleles: {a}, count: {b}, prc: {prc:.1f}%')