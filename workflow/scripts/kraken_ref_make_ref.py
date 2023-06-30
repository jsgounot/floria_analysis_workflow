# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-04-24 10:46:51
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-04-25 10:45:00

import gzip, re
from Bio import SeqIO

# Checking that the same record id is not used twice
RECORDIDS = set()

def slugify(value):
    # https://stackoverflow.com/questions/295135/turn-a-string-into-a-valid-filename
    # Slighlty MODIFIED from https://github.com/django/django/blob/master/django/utils/text.py
    value = re.sub(r'[^\w\s-]', '_', value.lower())
    return re.sub(r'[-\s]+', '-', value).strip('-_')

def sanitize_records(records):
    # here I sanitize record names for future usage
    # since record names will sometime end up into filepath as
    # folder name, it leads to annoying errors.
    for record in records:
        record.id = slugify(record.id)

        if record.id in RECORDIDS:
            raise Exception(f'Record id {record.id} already found once')

        RECORDIDS.add(record.id)
        yield record

def extract_seq(fname, tids):
    with gzip.open(fname, 'rt') as f:
        fdata = SeqIO.parse(f, 'fasta')
        for record in fdata:
            seqname, _, taxid = record.id.split('|')
            if int(taxid) in tids:
                yield record                

taxlist = snakemake.input['taxlist']
with open(taxlist) as f:
    taxlist = {int(line.strip()) for line in f}

fasta = snakemake.input['fasta']
outfile = snakemake.output[0]
with gzip.open(outfile, 'wt') as f:
    records = extract_seq(fasta, taxlist)
    records = sanitize_records(records)
    SeqIO.write(records, f, 'fasta')