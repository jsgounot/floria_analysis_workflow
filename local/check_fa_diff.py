from Bio import SeqIO

ref = '../res/ecoli_sim_basic/references/used/simulated_NC_000913_4/sample/ref.fasta'
ref = {record.id: record.seq for record in SeqIO.parse(ref, 'fasta')}

que = '../res/ecoli_sim_basic/references/used/simulated_NC_000913_4/sample/mut.fasta'
que = {record.id: record.seq for record in SeqIO.parse(que, 'fasta')}

diff = length = 0

for record, rseq in ref.items():
    for position, rbase in enumerate(rseq):
        if rbase != que[record][position]:
            diff += 1
        length += 1

print (diff, length, diff * 100 / length)

from collections import Counter

c = Counter(ref['NC_000913.3'])
print (c, sum(c.values()))