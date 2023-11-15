import os

fnames = snakemake.input
outfile = snakemake.output[0]
tempfile = outfile + '.temp'

print (len(fnames))

for i in range(0, len(fnames), 1000):
	sub = fnames[i:i+1000]
	sub = ' '.join(sub)
	if i == 0:
		cmdline = f'cat {sub} > {outfile}'
		os.system(cmdline)
	else:
		cmdline = f'cat {sub} > {tempfile}'
		os.system(cmdline)

		cmdline = f'cat {outfile} {tempfile} > {outfile}'
		os.system(cmdline)

		os.remove(tempfile)