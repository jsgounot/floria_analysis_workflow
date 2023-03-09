# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-03-25 12:02:50
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-03-25 12:09:38

import pysam
import sys

fname = sys.argv[1]
outfi = sys.argv[2]

fname = pysam.AlignmentFile(fname, 'rb')
outfi = pysam.AlignmentFile(outfi, 'wb', template=fname)

reads_ids = set()

for read in fname.fetch():
	if read.is_secondary or read.is_supplementary or read.query_name in reads_ids:
		continue

	reads_ids.add(read.query_name)
	outfi.write(read)

fname.close()
outfi.close()