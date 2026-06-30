import json
import pandas as pd

# Load R20 Top 6
with open(r'D:\workspace\round20\r20_top6.json') as f:
    top6 = json.load(f)

# Load Exclusion list
with open(r'D:\生信\2026Protein Design\Exclusion_List.csv') as f:
    excl = set(line.strip() for line in f if line.strip() and not line.startswith(','))

# Load R19 for compare
with open(r'D:\workspace\round19\final_6_r19.json') as f:
    r19 = json.load(f)

print('='*85)
print(' R20 TOP 6 ANALYSIS '.center(85, '='))
print('='*85)
print(f'{"Seq":>4s} {"Score":>7s} {"pTM":>7s} {"pLDDT":>7s} {"Chromo":>7s} {"Parent"}')
print('-'*85)
all_pass = True
for i, c in enumerate(top6):
    s = c['seq']
    m = s[0]=='M'
    L = 220 <= len(s) <= 250
    aa = all(a in 'ACDEFGHIKLMNPQRSTVWY' for a in s)
    excl_ok = s not in excl
    ok = m and L and aa and excl_ok
    if not ok: all_pass = False
    pn = c['parent'][:40]
    print(f'{i+1:>4d} {c["score"]:>7.4f} {c["ptm"]:>7.4f} {c["plddt"]:>7.3f} {c["chromo"]:>7.3f} {pn}')

print('-'*85)
print(f'Overall: {"✅ ALL PASS" if all_pass else "❌ SOME FAILED"}')
print()
print(f'R19 Top 1: {r19[0]["score"]:.4f}')
print(f'R20 Top 1: {top6[0]["score"]:.4f} (NEW RECORD!)')
delta = top6[0]['score'] - r19[0]['score']
print(f'Delta: {delta:+.4f} ({delta/r19[0]["score"]*100:+.2f}%)')

scores = [c['score'] for c in top6]
print(f'\nTop 6 range: {min(scores):.4f} ~ {max(scores):.4f}')
print(f'Avg: {sum(scores)/6:.4f}')
avg_r20 = sum(scores)/6
avg_r19 = sum(c['score'] for c in r19)/6
print(f'Avg improvement: {avg_r20 - avg_r19:+.4f}')

print(f'\nChromo pLDDT:')
print(f'  R19: {[c["chromo"] for c in r19]}')
print(f'  R20: {[c["chromo"] for c in top6]}')
