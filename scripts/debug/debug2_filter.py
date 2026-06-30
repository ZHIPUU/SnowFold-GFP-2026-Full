import json, os

# Check what MPNN produced vs parent
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

# parent is the FIRST entry (not the candidate)
parent_name, parent_seq = mpnn[0]
print(f'Parent from FA: {parent_name[:60]}')
print(f'Parent seq: {parent_seq[:60]}')
print(f'Parent seq len: {len(parent_seq)}')

# R19[0] parent seq
r = json.load(open('/root/autodl-tmp/r19/final_6_r19.json'))
ps = r[0]['seq']
print(f'\nR19[0] seq: {ps[:60]}')
print(f'R19[0] len: {len(ps)}')

print(f'\nMatch: {parent_seq == ps}')
if parent_seq != ps:
    # Find first difference
    for i in range(min(len(parent_seq), len(ps))):
        if parent_seq[i] != ps[i]:
            print(f'First diff at pos {i}: parent_fa={parent_seq[i]}, r19={ps[i]}')
            print(f'  context: ...{parent_seq[max(0,i-5):i+5]}... vs ...{ps[max(0,i-5):i+5]}...')
            break

# Check fixed positions (0-indexed = 64, 65, 66, 95, 221)
print('\nFixed positions (0-indexed): 64, 65, 66, 95, 221')
for pos in [64, 65, 66, 95, 221]:
    print(f'  pos {pos}: parent={parent_seq[pos]}, r19={ps[pos]}')

# Check candidates
print('\nFirst 5 candidates:')
for n, s in mpnn[1:6]:
    diff_count = sum(1 for a,b in zip(parent_seq, s) if a != b)
    print(f'  {n[:50]:50s} diff={diff_count:3d} first3={s[:3]}')

# Check: do candidates preserve fixed positions?
print('\nFixed position preservation check:')
fixed_preserved = 0
for n, s in mpnn[1:]:
    if all(s[p] == parent_seq[p] for p in [64, 65, 66, 95, 221]):
        fixed_preserved += 1
print(f'  {fixed_preserved}/{len(mpnn)-1} preserve all 5 fixed positions')

# Check: do candidates start with M?
m_count = sum(1 for n, s in mpnn[1:] if s.startswith('M'))
print(f'  {m_count}/{len(mpnn)-1} start with M')
