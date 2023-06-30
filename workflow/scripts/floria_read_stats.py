import tarfile
import fnmatch
import os
import gzip
import pandas as pd


def read_fastq(fdata):
    readid = None
    for idx, line in enumerate(fdata):
        if readid:
            line = line.decode("utf-8").strip()
            yield {'readid': readid, 'length': len(line)}
            readid = None
        if idx % 4 == 0:
            readid = line.decode("utf-8").split()[0].strip()

def iter_reads_info_tar(fname):
    pattern = '*.fastq'

    with tarfile.open(fname, "r:*") as tar:
        for fname in tar.getnames():
            bname = os.path.basename(fname)
            if fnmatch.fnmatch(fname, pattern):
                for readinfo in read_fastq(tar.extractfile(fname)):
                    yield {'bname': bname,  ** readinfo}


def iter_reads_info_gzips(fnames):
    for fname in fnames:
        bname = os.path.basename(fname)
        with gzip.open(fname) as f:
            for readinfo in read_fastq(f):
                yield {'bname': bname,  ** readinfo}

res = snakemake.input['res']
res = pd.DataFrame(iter_reads_info_tar(res))
res['cat'] = 'phased'

allreads = snakemake.input['reads']
allreads = pd.DataFrame(iter_reads_info_gzips(allreads))
allreads = allreads[~ allreads['readid'].isin(res['readid'])]
allreads['cat'] = 'unphased'

df = pd.concat([res, allreads])
outfile = snakemake.output[0]
df.to_csv(outfile, sep='\t')

#df['sample'] =  df['readid'].apply(lambda rid: rid.split('-')[0][1:])
#df = df.groupby(['bname', 'sample'])['length'].agg(['size', 'sum'])
#print (df)

'''
print (res)
print (allreads)
print (df)
print (df['cat'].value_counts())
print (res.iloc[0].to_dict())
print (allreads.iloc[0].to_dict())
'''