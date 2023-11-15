import tarfile, gzip

import logging
from logging.handlers import RotatingFileHandler

from pybedtools import BedTool, Interval
from pyfaidx import Fasta, FastaRecord
import pysam

# Logger
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s')
 
logpath = snakemake.log[0]
file_handler = RotatingFileHandler(logpath, 'w')
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

# stream_handler = logging.StreamHandler()
# stream_handler.setFormatter(formatter)
# stream_handler.setLevel(logging.DEBUG)
# logger.addHandler(stream_handler)

fasta = snakemake.input['fasta']
vcf = snakemake.input['vcf']
archive = snakemake.input['reads_ar']
outfile = snakemake.output[0]
default_ref = bool(snakemake.params['default_ref'])
add_nocov = bool(snakemake.params['add_nocov'])

# fasta =   'res/synthetic_com_2_2/temp/test/132.illumina.fa'
# vcf =     'res/synthetic_com_2_2/temp/test/132.illumina.freebayes.floria_vcf_header.new.vcf.gz'
# archive = 'res/synthetic_com_2_2/temp/test/res.tar.gz'
# outfile = 'res/synthetic_com_2_2/temp/test/out.new.fa.gz'
# default_ref = True
# add_nocov = True

def parse_vartig(stream):
    r = []
    for line in stream:
        hapid, contig_name, snprange, baserange, cov, err, hapq, rel_err = line.decode('utf-8').split('\t')
        contig_name = contig_name.split(':')[1]
        snprange = tuple(map(int, snprange.split(':')[1].split('-')))
        baserange = tuple(map(int, baserange.split(':')[1].split('-')))
        vartig = next(stream).decode('utf-8').strip()
        r.append((hapid, contig_name, snprange, baserange, vartig))

    return r

d = {}
tar = tarfile.open(archive)
for member in tar.getmembers():
    if member.name.endswith('.vartigs'):
        logger.info(f'Parse: {member.name}')
        f = tar.extractfile(member)
        d[member.name] = parse_vartig(f)


logger.info(f'Parse: {fasta}')
fasta = Fasta(fasta)

logger.info(f'Parse: {vcf}')
vcf = pysam.VariantFile(vcf)

with gzip.open(outfile, 'wt') as of:
    for contig, cinfo in d.items():
        for haplotig in cinfo:
            hapid, contig_name, snprange, baserange, vartig = haplotig
            logger.debug(f'Process: {haplotig}')

            sequence = fasta.get_seq(contig_name, baserange[0], baserange[1]).seq
            sequence = list(sequence)  # convert string to list for mutability

            startpos = baserange[0]
            variants = vcf.fetch(contig_name, startpos - 1, baserange[1])
            variants = {variant.pos: [variant.ref] + list(variant.alts) 
                for variant in variants 
                if len(variant.ref) == 1 and all(len(alt) == 1 for alt in variant.alts)
                }

            logger.debug(f'{len(variants)}, {len(vartig)}, {baserange}, {vartig}')
            assert len(variants) == len(vartig)

            for idx, position in enumerate(sorted(variants)):
                seqposition = position - startpos
                alts = variants[position] # alts[0] == ref
                
                # Asserting potential issues
                if sequence[seqposition] != alts[0]:
                    ref = sequence[seqposition].upper()
                    ref = ref if ref in 'ATGC' else 'N'
                    if ref != alts[0]:
                        emsg = f'{seqposition}, {idx}, {position}, {sequence[seqposition]} {alts[0]}, {alts}'
                        raise Exception(emsg) 

                vartigvalue = vartig[idx]
                # vartigvalue = 0 for ref, ? for not enough covered, or 1,2,... for alts

                if vartigvalue == '0': continue
                elif vartigvalue == '?':
                    if default_ref: continue
                    else: sequence[seqposition] = 'N'
                else: sequence[seqposition] = alts[int(vartigvalue)]

            sequence = ''.join(sequence)
            hapid = contig_name + '_' + hapid.split('.')[0][1:]
            of.write(f'>{hapid}\n')
            of.write('\n'.join(sequence[i:i+60] for i in range(0, len(sequence), 60)) + '\n')

    if add_nocov:
        # Reference genome parts not covered by vartigs
        frange = BedTool((Interval(record.name, 0, len(record)) for record in fasta))
        vrange = BedTool((Interval(haplotig[1], haplotig[3][0], haplotig[3][1], name=haplotig[0]) 
            for haplogroup in d.values() for haplotig in haplogroup)).merge()
        subtract = frange.subtract(vrange)

        for idx, interval in enumerate(subtract):
            hapid = f'>{interval.chrom}_nocov_{idx}'
            sequence = fasta.get_seq(interval.chrom, 
                interval.start + 1, interval.end).seq
            
            of.write(f'>{hapid}\n')
            of.write('\n'.join(sequence[i:i+60] for i in range(0, len(sequence), 60)) + '\n')