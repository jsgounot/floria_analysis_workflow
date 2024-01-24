# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-08-10 14:23:06
# @Last Modified by:   jsgounot
# @Last Modified time: 2024-01-23 09:50:48

import json

NAME_ATTRIBUTES = {
    'kraken_ref': ('mode', ),
    'kraken_ref_merged': ('mode', ),
    'flye': ('read', 'mode'),
    'floria': ('readtype', 'fmode', 'post_assembler', 'assembler_rtype', 'assembler_preset'),
    'floria_single': ('readtype', 'fmode'),
    'strainxpress': ('mode',),
    'strainberry': ('readtype',),
    'megahit': (),
    'inpref': ()
}

PHASER_REQUIREMENTS = {
    'floria': ('vcalling', 'assembly'),
    'floria_single': ('vcalling', 'assembly'),
    'strainberry': ('assembly', )
}

UNDERSCORE_SEP = (
    'flye',
    )

DIRECT_READS_PHASERS = (
    'strainxpress',
    )

class ConfigParserError(Exception):
    pass

def get_raise(dic, key):
    try: 
        return dic[key]
    except KeyError:
        msg = f'Unable to find mandatory "{key}" key in json file, please refer to the documentation.'
        raise ConfigParserError(msg)

def parse_sub(param, mainkey, optional=False, ignore_attributes=(), sep='.'):
    try: 
        subparam = get_raise(param, mainkey)
    except ConfigParserError as e:
        if optional:
            return mainkey, ''
        else:
            raise e

    name = subpath = get_raise(subparam, 'name')

    # TODO: This need to be removed, but will need to change snakemake pathnames ...
    sep = '_' if name in UNDERSCORE_SEP else '.'

    attributes = NAME_ATTRIBUTES.get(name, None)
    if attributes is None: 
        raise ConfigParserError(f'Unable to find attributes linked to name "{name}", please refer to the documentation.')

    for attribute in attributes:
        if attribute in ignore_attributes: continue
        value = get_raise(subparam, attribute)
        subpath += sep + value

    return name, subpath

def parse_rfiltering(param):
    attributes = ('similarity', 'length')
    subparam = get_raise(param, 'ref_filtering')

    sim = get_raise(subparam, 'similarity')
    if sim != 'none': 
        try: 
            float(sim)
            sim = sim.replace('.', '')
        except ValueError:
            ConfigParserError(f'Similarity value in ref_filtering must be a number or "none", not "{sim}"')

    length = get_raise(subparam, 'length')
    if length != 'none':
        try :
            int (length)
        except ValueError:
            ConfigParserError(f'Length value in ref_filtering must be an integer or "none", not "{length}"')

    return f'done.{sim}.{length}.empty'

def parse_refcomp(values, groups):
    res = []
    for value in values:
        value = value.lower()
        if value == 'mummer':
            res += [f'refcomp/{group}.mummer.txt' for group in groups]
        elif value == 'fastani':
            res += [f'refcomp/{group}.fastani.txt' for group in groups]
        else:
            raise ConfigParserError(f'refcomp values must be either "mummer" or "fastani", found "{value}"')

    return res

def parse_config(jdata, groups):
    outputs = []

    for param in jdata:
        try:
            aname, asub = parse_sub(param, 'assembly')
        except ConfigParserError as e:
            pname, psub = parse_sub(param, 'phasing')
            if pname in DIRECT_READS_PHASERS:
                aname = asub = ''
            else:
                raise e

        if aname == 'kraken_ref':
            pname, psub = parse_sub(param, 'phasing', True, ('readtype',))
        else:
            pname, psub = parse_sub(param, 'phasing', True)

        if pname in ('floria', ):
            vsub = get_raise(param, 'vcalling')
        else:
            vsub = ''

        phasename = asub
        if vsub: phasename += '.' + vsub
        if psub: phasename += '.' + psub

        basename = parse_rfiltering(param)

        group = get_raise(param, 'group')
        used_groups = groups if group == 'all' else [group]
        for group in used_groups:
            fname = f'stats/assemblies/{group}/{phasename}/mummer/circos/{basename}'
            outputs.append(fname)

        refcomp = param.get('refcomp', [])
        outputs.extend(parse_refcomp(refcomp, used_groups))

    return outputs

def parse_config_production(jdata, groups):
    outputs = []

    for param in jdata:
        try:
            aname, asub = parse_sub(param, 'assembly')
        except ConfigParserError as e:
            pname, psub = parse_sub(param, 'phasing')
            if pname in DIRECT_READS_PHASERS:
                aname = asub = ''
            else:
                raise e

        if aname == 'kraken_ref':
            pname, psub = parse_sub(param, 'phasing', True, ('readtype',))
            pname, ssub = parse_sub(param, 'phasing', True, ('readtype', 'post_assembler', 'assembler_rtype', 'assembler_preset'))
        else:
            pname, psub = parse_sub(param, 'phasing', True)
            pname, ssub = parse_sub(param, 'phasing', True, ('post_assembler', 'assembler_rtype', 'assembler_preset'))

        if pname in ('floria', 'floria_single'):
            vsub = get_raise(param, 'vcalling')
        else:
            vsub = ''

        make_basename = lambda l, sep: sep.join(e for e in l if e)

        if pname == 'floria_single':
            psub = psub.replace('floria_single', 'floria')
            ssub = ssub.replace('floria_single', 'floria')

        group = get_raise(param, 'group')
        used_groups = groups if group == 'all' else [group]
        for group in used_groups:

            if pname in ('floria', 'floria_single'):
                basename = make_basename((asub, vsub, ssub), '.')
                fname = f'results/{group}/phasing/floria/{basename}.tsv.gz'
                outputs.append(fname)
                
            if pname in ('floria', 'strainberry'):
                basename = make_basename((asub, vsub, psub), '.')
                fname = f'results/{group}/phasing/{pname}/{basename}.fa.gz'
                outputs.append(fname)

        refcomp = param.get('refcomp', [])
        outputs.extend(parse_refcomp(refcomp, used_groups))

    return outputs

'''
phasing/floria/split/{group}/{krtype}/{tid}/{vcaller}/ploidy.{fmode}.tsv.gz
expand('phasing/floria/presplit/{group}/nanopore/longshot/ploidy.tsv.gz', group=GROUPS),
expand('phasing/floria/presplit/{group}/nanopore/longshot/assembly.wtdbg2.long_reads.nano.fa.gz')
'temp/assemblies/{group}/kraken_presplit.{krtype}.{vcaller}.floria.{assembler}.{rtype}.{preset}.fa.gz'
'''