import concurrent.futures
import glob, sys, os
import pandas as pd
import pysam

import depthcomp

dname = sys.argv[1]

bamfiles = os.path.join(dname, 'mapping/*/*.bam')
bamfiles = glob.glob(bamfiles)
assert bamfiles

def process_bam(bamfile):
	print (f'Read: {bamfile}')
	df = depthcomp.make_depth_window(bamfile)
	return {
		'fname': bamfile,
		'npos': df['npos'].sum(),
		'mcov': df['sum_depth'].sum() / df['npos'].sum()
	}

df = []
with concurrent.futures.ProcessPoolExecutor() as executor:
	for res in executor.map(process_bam, bamfiles):
		df.append(res)

bamdata = pd.DataFrame(df)

def process_vcf(fname):
	print (f'Read: {fname}')
	vcf = pysam.VariantFile(fname)
	return {
		'fname': fname,
		'snpcount': sum(1 if rec.ref != rec.alts else 0 for rec in vcf.fetch())
	}

vcffiles = os.path.join(dname, 'vcalling/*/*.vcf')
vcffiles = glob.glob(vcffiles)
assert vcffiles

df = []
with concurrent.futures.ProcessPoolExecutor() as executor:
	for res in executor.map(process_vcf, vcffiles):
		df.append(res)

vcfdata = pd.DataFrame(df)


fun = lambda fname: os.path.basename(os.path.dirname(fname))
bamdata['sample'] = bamdata['fname'].apply(fun)
vcfdata['sample'] = vcfdata['fname'].apply(fun)

fun = lambda fname: os.path.splitext(os.path.basename(fname))[0]
bamdata['bname'] = bamdata['fname'].apply(fun)
vcfdata['bname'] = vcfdata['fname'].apply(fun)
vcfdata['bname'] = vcfdata['bname'].apply(fun)

df = bamdata.merge(vcfdata, on=['sample', 'bname'], how='outer')
outfile = os.path.join(dname, 'reports/mappingstat.tsv.gz')
df.to_csv(outfile, sep='\t', index=False, compression='gzip')