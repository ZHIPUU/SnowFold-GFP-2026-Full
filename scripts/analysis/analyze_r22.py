import json, csv

# Load R22 Top 6
with open(r'D:\workspace\round22\final_6_r22.json') as f:
    top6 = json.load(f)

# Load Exclusion list
with open(r'D:\生信\2026Protein Design\Exclusion_List.csv') as f:
    excl = set(line.strip() for line in f if line.strip() and not line.startswith(','))

# Compare with R19
with open(r'D:\workspace\round19\final_6_r19.json') as f:
    r19 = json.load(f)

# Compare with R20
with open(r'D:\workspace\round20\r20_top6.json') as f:
    r20 = json.load(f)

print('='*100)
print(' R22 TOP 6 ANALYSIS (NEW PROJECT RECORD) '.center(100, '='))
print('='*100)
print(f"{'Seq':>4s} {'Score':>7s} {'pTM':>7s} {'pLDDT':>7s} {'Chromo':>7s} {'r':>3s} {'Parent'}")
print('-'*100)
all_pass = True
for i, c in enumerate(top6):
    s = c['seq']
    m_ok = s[0] == 'M'
    len_ok = 220 <= len(s) <= 250
    aa_ok = all(a in 'ACDEFGHIKLMNPQRSTVWY' for a in s)
    excl_ok = s not in excl
    ok = m_ok and len_ok and aa_ok and excl_ok
    if not ok: all_pass = False
    pn = c.get('parent', 'unknown')[:35]
    r = c.get('recycles', 8)
    print(f"{i+1:>4d} {c['score']:>7.4f} {c['ptm']:>7.4f} {c['plddt']:>7.3f} {c['chromo']:>7.3f} {r:>3d} {pn}")

print('-'*100)
print(f'Overall: {"ALL PASS" if all_pass else "SOME FAILED"}')
print()
print(f'R19 Top 1: {r19[0]["score"]:.4f}')
print(f'R20 Top 1: {r20[0]["score"]:.4f}')
print(f'R22 Top 1: {top6[0]["score"]:.4f} (NEW RECORD)')
delta_r19 = top6[0]['score'] - r19[0]['score']
delta_r20 = top6[0]['score'] - r20[0]['score']
print(f'Delta vs R19: {delta_r19:+.4f} ({delta_r19/r19[0]["score"]*100:+.2f}%)')
print(f'Delta vs R20: {delta_r20:+.4f} ({delta_r20/r20[0]["score"]*100:+.2f}%)')

# Avg
avg_r22 = sum(c['score'] for c in top6) / 6
avg_r20 = sum(c['score'] for c in r20) / 6
avg_r19 = sum(c['score'] for c in r19) / 6
print()
print(f'Avg Top 6:')
print(f'  R19: {avg_r19:.4f}')
print(f'  R20: {avg_r20:.4f} ({avg_r20 - avg_r19:+.4f})')
print(f'  R22: {avg_r22:.4f} ({avg_r22 - avg_r20:+.4f})')
