# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-04-02 23:46:46
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-04-02 23:47:56

import pandas as pd
import glob, os

fnames = 'stats/assemblies/*/*/mummer/report.tsv'
fnames = glob.glob(fnames)

for fname in fnames:
	print (fname)
	print (pd.read_csv(fname, sep='\t').drop('refpath', axis=1))