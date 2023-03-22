import sys, subprocess, io
import pandas as pd

def process_chunk(sdf, windowsize):
    sdf['window'] = (sdf['position'] - 1) // windowsize
    sdf = sdf.groupby(['contig', 'window'])['depth'].agg(['size', 'sum']).reset_index()
    sdf.columns = ['contig', 'window', 'npos', 'sum_depth']
    return sdf.reset_index()

def make_depth_window(bamfile, windowsize=10000):
    cmdline = f'samtools depth -aa {bamfile}'
    sp = subprocess.Popen(cmdline.split(), stdout=subprocess.PIPE)
    sio = io.TextIOWrapper(sp.stdout)
    names = ['contig', 'position', 'depth']
    iterator = pd.read_csv(sio, sep='\t', names=names, iterator=True, chunksize=100000)
    df = pd.concat(process_chunk(sdf, windowsize) for sdf in iterator)
    df = df.groupby(['contig', 'window'])['npos', 'sum_depth'].sum()
    df['mcov'] = df['sum_depth'] / df['npos']
    return df.reset_index()

if __name__ == '__main__':
    bamfile = sys.argv[1]
    try: windowsize = int(sys.argv[2])
    except IndexError: windowsize = 10000
    df = make_depth_window(bamfile, windowsize)
    print (df)  