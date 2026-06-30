import json

with open(r'D:\workspace\round23\final_6_r23.json') as f: r23 = json.load(f)
with open(r'D:\workspace\round22\final_6_r22.json') as f: r22 = json.load(f)
with open(r'D:\workspace\round20\r20_top6.json') as f: r20 = json.load(f)
with open(r'D:\workspace\round19\final_6_r19.json') as f: r19 = json.load(f)

print('='*80)
print('R19 -> R20 -> R22 -> R23 EVOLUTION')
print('='*80)
for nm, t in [('R19',r19),('R20',r20),('R22',r22),('R23',r23)]:
    s = t[0]["score"]; p = t[0]["ptm"]; l = t[0]["plddt"]; c = t[0]["chromo"]
    print(f'{nm} Top 1: {s:.4f} (pTM={p:.4f}, pLDDT={l:.3f}, chromo={c:.3f})')

print()
print('R23 Top 6:')
for i, c in enumerate(r23):
    pn = c["parent"]
    print(f'  {i+1}. score={c["score"]:.4f} pTM={c["ptm"]:.4f} parent={pn}')

print()
print('R22 Top 6 (current record):')
for i, c in enumerate(r22):
    print(f'  {i+1}. score={c["score"]:.4f} pTM={c["ptm"]:.4f}')

print()
# R23 Top 1 vs R22 Top 1
print('=== R23 vs R22 Top 1 (head-to-head) ===')
r23_top1 = r23[0]['score']
r22_top1 = r22[0]['score']
print(f'R23 Top 1: {r23_top1:.4f}')
print(f'R22 Top 1: {r22_top1:.4f} (RECORD)')
print(f'Delta: {r23_top1 - r22_top1:+.4f} ({(r23_top1-r22_top1)/r22_top1*100:+.2f}%)')

# Are R23 Top 6 different from R22 Top 6?
print()
print('=== Sequence overlap (R23 vs R22) ===')
r23_seqs = set(c['seq'] for c in r23)
r22_seqs = set(c['seq'] for c in r22)
overlap = r23_seqs & r22_seqs
print(f'R23 unique seqs: {len(r23_seqs)}')
print(f'R22 unique seqs: {len(r22_seqs)}')
print(f'Overlap: {len(overlap)}')
for s in overlap:
    print(f'  Common: {s[:50]}...')

# Best combined top 6 (R22 + R23 deduplicated)
print()
print('=== Combined Top 6 (R22 + R23, dedup) ===')
combined = r22 + r23
seen = set(); unique = []
for c in combined:
    if c['seq'] not in seen:
        seen.add(c['seq']); unique.append(c)
unique.sort(key=lambda x: x['score'], reverse=True)
for i, c in enumerate(unique[:6]):
    print(f'  {i+1}. score={c["score"]:.4f} from={c["parent"]}')
