# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-04-01 17:34:35
# @Last Modified by:   jsgounot
# @Last Modified time: 2024-03-14 14:39:11

DEFAULT_GENERIC_THREADS = 32
DEFAULT_GENERIC_THREADS_SPLIT = 16

import os, copy
import gzip

def get_ref(wildcards, config, scon=None):
    group, sample = wildcards.group, wildcards.sample   
    scon = scon or config['samples'][group][sample]

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
    scon = config['samples'][group][sample]

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
    sample = [sample for sample, sdata in config['samples'][group].items() 
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

def get_generic_threads(wc):
    return GENERIC_THREADS

def get_generic_threads_split(wc):
    return GENERIC_THREADS_SPLIT

def get_reads(group, mtype, model=None):
    if mtype == 'illumina':
        return [
            f'reads/{group}/illumina/merged/R1.fastq.gz',
            f'reads/{group}/illumina/merged/R2.fastq.gz'
        ]

    elif mtype == 'pacbio':
        return [f'reads/{group}/pacbio/merged_{model}.fastq.gz']

    elif mtype == 'nanopore' or mtype == 'hybrid':
        return [f'reads/{group}/nanopore/merged.fastq.gz']
    
    else:
        raise Exception(f'mtype not found: {mtype}')

def get_softpath(config, softname):
    softpath = config['softparams']['soft'].get(softname, None)
    if softpath is None:
        raise Exception(f'Configuration error: Unable to find software path for soft "{softname}" within softpaths.json')
    if softpath == '':
        print (f'WARNING: Softpath for soft name {softname} is currently empty. Please add it to your softpaths.json if you plan to use it.')
    return softpath

def get_vcalling_ploidy(config, wc, ploidy=None):
    ploidy = ploidy or wc.ploidy
    if ploidy == 'auto':
        ploidy = len(config['samples'][wc.group])

    try:
        int(ploidy)
    except ValueError:
        raise Exception(f'Ploidy "{ploidy}" cannot be converted to integer')

    if int(ploidy) > 5:
        print (f'Warning, ploidy for variants calling is higher than 5 ("{ploidy}"), are you sure?')

    return ploidy

def infer_vcalling_ploidy(config, wc):
    # infer vcalling ploidy from sample group name
    ploidy = wc.vcaller.split('_')[-1]
    return get_vcalling_ploidy(config, wc, ploidy)