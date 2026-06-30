import json, os
r = json.load(open('/root/autodl-tmp/r19/final_6_r19.json'))
print('parent[0] name:', r[0]['name'][:50])
print('parent[0] seq:', r[0]['seq'][:60])
print('len:', len(r[0]['seq']))

mpnn = []
fa = '/root/autodl-tmp/r20/mpnn_out/T=0.2, sample=25, score=0.6229, global_score=0.6383, seq_recovery=0.8291/seqs/T=0.2, sample=25, score=0.6229, global_score=0.6383, seq_recovery=0.8291.fa'
with open(fa) as f:
    name=''; seq=''
    for line in f:
        line=line.strip()
        if line.startswith('>'):
            if seq: mpnn.append((name, seq))
            name=line[1:]; seq=''
        else:
            if line: seq+=line
    if seq: mpnn.append((name, seq))
print('total:', len(mpnn))
ps = r[0]['seq']
match = sum(1 for n,s in mpnn if s==ps)
nonM = sum(1 for n,s in mpnn if not s.startswith('M'))
unique_len = len(set(s for n,s in mpnn))
print(f'==parent: {match}, nonM: {nonM}, unique: {unique_len}')
filtered = [s for n,s in mpnn if s!=ps and s.startswith('M')]
print('filtered (valid):', len(filtered))
