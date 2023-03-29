# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-01-19 15:35:44
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-01-19 15:44:03

import glob, os, sys
from itertools import combinations
from Bio import SeqIO

dname = sys.argv[1]
groups = os.path.join(dname, 'references/used/*/')
groups = sorted(glob.glob(groups))

for group in groups:
	print (f'{os.path.basename(os.path.dirname(group))}:')
	fnames = os.path.join(group, 'sample/*.fasta')
	fnames = glob.glob(fnames)
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
					tot += 1

		if tot:
			prc = diff * 100 / tot
			print (f'{os.path.basename(a)} {os.path.basename(b)} {diff} {tot} {prc:.1f}')
		else:
			print (f'{os.path.basename(a)} {os.path.basename(b)} Nothing to compare')