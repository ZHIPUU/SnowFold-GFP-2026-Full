#!/usr/bin/env python3
"""Check R19 submission compliance"""
import json

# Load R19 final 6
with open(r'D:\workspace\round19\final_6_r19.json') as f:
    top6 = json.load(f)

# Load exclusion list
with open(r'D:\生信\2026Protein Design\Exclusion_List.csv') as f:
    excl = set(line.strip() for line in f if line.strip() and not line.startswith(','))

print("="*80)
print(f"{'#':<4s} {'Name':<58s} {'Score':<8s} {'M?':<4s} {'Len':<5s} {'AA?':<4s} {'Excl?':<6s} {'Result':<8s}")
print("="*80)

all_pass = True
for i, c in enumerate(top6):
    s = c['seq']
    m_ok = s[0] == 'M'
    len_ok = 220 <= len(s) <= 250
    aa_ok = all(a in 'ACDEFGHIKLMNPQRSTVWY' for a in s)
    excl_ok = s not in excl
    passed = m_ok and len_ok and aa_ok and excl_ok
    if not passed: all_pass = False
    name = c['name'][:55]
    print(f"{i+1:<4d} {name:<58s} {c['score']:<8.4f} {'✓' if m_ok else '✗':<4s} {len(s):<5d} {'✓' if aa_ok else '✗':<4s} {'✓' if excl_ok else '✗':<6s} {'✅' if passed else '❌':<8s}")

print("="*80)
print(f"\nOverall: {'ALL PASS ✅' if all_pass else 'SOME FAILED ❌'}")
print(f"Scores: {[c['score'] for c in top6]}")
print(f"Range: {min(c['score'] for c in top6):.4f} ~ {max(c['score'] for c in top6):.4f}")
