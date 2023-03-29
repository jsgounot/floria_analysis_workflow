# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-03-28 14:44:03
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-03-28 15:18:23

import json

gradient = list(range(81, 96, 2))
jdata = {}

for idx, error_rate in enumerate(gradient):

	name = f'error_rate_{idx}'

	jdata[name] = {
		'ASM1990034v1': {
			"ncbiasbly": 'ASM1990034v1',
			"mapping_ref": "True",
			"quantity": 15,
			"circular": "true",
			"nanopore": f"average_pid{error_rate}",
			"seed": 1
		},
		'ASM1990032v1': {
			"ncbiasbly": 'ASM1990032v1',
			"quantity": 30,
			"circular": "true",
			"nanopore": f"average_pid{error_rate}",
			"seed": 1
		},
		'ASM1990030v1': {
			"ncbiasbly": 'ASM1990030v1',
			"quantity": 20,
			"circular": "true",
			"nanopore": f"average_pid{error_rate}",
			"seed": 1
		}
	}

outfile = './kleb_error_rate.json'
with open(outfile, 'w') as f:
	json.dump(jdata, f, indent=4)