from Bio import SeqIO
import os

# Need to add the hidx so contig names are not the same
# when assembly are merged together

fnames = snakemake.input
outfile = snakemake.output[0]

def iter_records(fname):
	hidx = os.path.basename(os.path.dirname(fname))
	hidx = hidx.split('.')[1].split('_')[1]
	# print (fname, hidx)
	for record in SeqIO.parse(fname, 'fasta'):
		record.id = record.id + '_' + hidx
		record.description = fname
		yield record

if not any(record for fname in fnames for record in iter_records(fname)):
	raise Exception('Unable to find any contigs. Either all haplotypes assembly failed or sample \
	    is completly unable to be phased. Note that some assemblers failures are ignored to prevent \
	    pipeline break by small haplotypes.')

records = (record for fname in fnames for record in iter_records(fname) if os.path.isfile(fname))
SeqIO.write(records, outfile, 'fasta')