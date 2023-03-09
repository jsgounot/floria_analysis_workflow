# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-01-19 15:35:44
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-01-19 15:44:03

import copy
import json

template  = '../configs/klebsiella.json'
sample    = 'klebsiella'
pid_range = range(80, 95, 1)
outfile   = '../configs/klebsiella_pid.json'

with open(template) as f:
	jdata = json.load(f)


sdata = jdata[sample]
ndata = {}

for pid in pid_range:
	sub = f'{sample}_pid{pid}'
	ndata[sub] = copy.deepcopy(sdata)

	for sample, samp_data in ndata[sub].items():
		samp_data['nanopore'] = f'average_pid{pid}'

with open(outfile, 'w') as f:
	json.dump(ndata, f, indent=4, sort_keys=True)