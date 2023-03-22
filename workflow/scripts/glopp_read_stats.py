import tarfile
import fnmatch
import os
import pandas as pd


def read_fastq(fdata):
    readid = None
    for idx, line in enumerate(fdata):
        if readid:
            line = line.decode("utf-8").strip()
            yield {'readid': readid, 'length': len(line)}
            readid = None
        if idx % 4 == 0:
            readid = line.decode("utf-8").strip()

def iter_reads_info(fname):
    pattern = '*.fastq'

    with tarfile.open(fname, "r:*") as tar:
        for fname in tar.getnames():
            bname = os.path.basename(fname)
            if fnmatch.fnmatch(fname, pattern):
                for readinfo in read_fastq(tar.extractfile(fname)):
                    yield {'bname': bname,  ** readinfo}


fname = snakemake.input[0]
df = pd.DataFrame(iter_reads_info(fname))
outfile = snakemake.output[0]
df.to_csv(outfile, sep='\t')

#df['sample'] =  df['readid'].apply(lambda rid: rid.split('-')[0][1:])
#df = df.groupby(['bname', 'sample'])['length'].agg(['size', 'sum'])
#print (df)