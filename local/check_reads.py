import gzip

s1, s2 = set(), set()
counter = 0

fname = '/mnt/volume1/strainphasing/glopp_analysis_workflow/res/kleb_covgrad_3strains/reads/simulated_kleb_cg_7_2/nanopore/merged.fastq.gz'
with gzip.open(fname, 'rb') as f:
	for idx, line in enumerate(f):
		if idx % 4 == 0:
			line = line.decode('utf-8').strip()
			counter += 1
			s1.add(line)
			s2.add(line.split()[0]) 

print (counter, len(s1), len(s2))