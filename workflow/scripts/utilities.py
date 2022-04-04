# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-04-01 17:34:35
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-04-01 17:36:12

def get_mapping_ref(wildcards, config):
    group = wildcards.group
    sample = [sample for sample, sdata in config[group].items() 
              if sdata.get('mapping_ref', None) == 'True']
    
    if len(sample) != 1:
        raise Exception(f'No or more than one sample is defined for reference mapping, group: {group}')

    sample = sample[0]

    return f'references/used/{group}/{sample}.fasta'