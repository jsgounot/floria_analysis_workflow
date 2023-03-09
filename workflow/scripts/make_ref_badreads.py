from Bio import SeqIO
import gzip

config = snakemake.params[0]
outfile = snakemake.output[0]
records = []

for sample, sdata in config.items():
	if sample == 'internal':
		continue

	fdata = SeqIO.parse(sdata['reference'], 'fasta')

	quantity = sdata['quantity']
	ratio = round(int(quantity) / 10, 1)

	circular = sdata.get('circular', None)
	circular = 'circular=' + circular if circular is not None else ''

	for idx, record in enumerate(fdata):
		record.id = f'{sample}_{idx} depth={ratio} {circular}'
		record.name = record.description = ''
		records.append(record)

with open(outfile, 'w') as f:
	SeqIO.write(records, f, 'fasta')