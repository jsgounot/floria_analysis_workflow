# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-04-01 17:34:35
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-03-06 10:26:47

import os, copy
import gzip

def get_ref(wildcards, config, scon=None):
    group, sample = wildcards.group, wildcards.sample   
    scon = scon or config[group][sample]

    if 'simulate' in scon:
        return f'references/simulated/{group}/{sample}.simseq.genome.fa'
    elif 'refpath' in scon:
        return os.path.relpath(scon['refpath'])
    elif 'ncbiasbly' in scon:
        ncbiid = scon['ncbiasbly']
        return f'references/download/{group}/ncbiasbly_{ncbiid}.fasta'
    elif 'ncbinuc' in scon:
        ncbiid = scon['ncbinuc']
        return f'references/download/{group}/ncbinuc_{ncbiid}.fasta'
    else:
        raise Exception(f'Neither a reference file or a NCBI ID is provided for sample {sample} in group {group}')

def get_sim_param(wildcards, config, param):
    group, sample = wildcards.group, wildcards.sample   
    scon = config[group][sample]

    if param == 'ref':
        scon = copy.deepcopy(scon)
        scon.pop('simulate')
        return get_ref(wildcards, config, scon)

    if param == 'snp_count':
        ref = get_sim_param(wildcards, config, 'ref')
        refsize = get_ref_size(ref)
        snp_prc = get_sim_param(wildcards, config, 'snp_prc')
        return int(refsize * (snp_prc / 100))

    try:
        return scon['simulate'][param]
    except KeyError:
        raise Exception(f'Not able to retrieve simulate param {param} for {sample}')

def get_ref_size(refpath):
    with open(refpath) as f:
        return sum(
            len(line.strip())
            for line in f
            if not line.startswith('>')
            )

def get_mapping_ref(wildcards, config):
    group = wildcards.group
    sample = [sample for sample, sdata in config[group].items() 
              if sdata.get('mapping_ref', None) == 'True']
    
    if len(sample) != 1:
        raise Exception(f'No or more than one sample is defined for reference mapping, group: {group}')

    sample = sample[0]

    return f'references/used/{group}/sample/{sample}.fasta'

def is_nonempty_gz_file(name):
    with gzip.open(name, 'rb') as f:
        try:
            file_content = f.read(1)
            return len(file_content) > 0
        except:
            return False
