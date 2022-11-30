# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-11-23 16:57:05
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-11-24 13:20:39

import json

outfile = snakemake.output[0]
config = snakemake.params.config

with open(outfile, 'w') as f:
	json.dump(config, f, indent=4, sort_keys=True)