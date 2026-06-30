"""Round 2: Build full sequences from GFP_data.xlsx mutation descriptions.

Strategy:
1. Load WT sequences for 4 GFP types
2. Parse mutation strings (e.g. "A109D:N145D:I187V")
3. Apply mutations to WT to reconstruct actual sequences
4. Save full sequence dataframe for downstream ESM embedding

Output: work/round2/all_variants.parquet (seq, type, brightness)
"""
import re
import time
import pandas as pd
import numpy as np
from pathlib import Path

WORK = Path(r"D:\生信\2026Protein Design")
ROUND2 = WORK / "work" / "round2"
ROUND2.mkdir(exist_ok=True, parents=True)

# ---- Load WT sequences ----
print("="*70)
print("Loading WT sequences...")
print("="*70)
WT = {}
with open(WORK / "AAseqs of 5 GFP proteins_20260511.txt") as f:
    cur_name = None
    cur_seq = []
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if line.startswith('>'):
            if cur_name and cur_seq:
                WT[cur_name] = ''.join(cur_seq)
            cur_name = line[1:].strip()
            cur_seq = []
        else:
            cur_seq.append(line)
    if cur_name and cur_seq:
        WT[cur_name] = ''.join(cur_seq)

print(f"Loaded {len(WT)} WT sequences:")
for n, s in WT.items():
    print(f"  {n}: {len(s)} aa")

# ---- Load GFP_data.xlsx ----
print("\nLoading GFP_data.xlsx...")
t0 = time.time()
df = pd.read_excel(WORK / "GFP_data.xlsx")
print(f"  {len(df):,} rows in {time.time()-t0:.1f}s")

# ---- Parse mutations ----
print("\nParsing mutations and applying to WT...")
t0 = time.time()
MUT_RE = re.compile(r'^([A-Z])(\d+)([A-Z])$')

def parse_mut(m):
    """Parse single mutation like 'A109D' -> (0-indexed pos, new_aa)."""
    m = MUT_RE.match(m.strip())
    if not m:
        return None
    return int(m.group(2)) - 1, m.group(3)  # 0-indexed

def apply_mutations(wt_seq, mut_str):
    """Apply mutations to WT sequence. Returns new seq or None on error.

    Accepts separators: '/', ':', or ','
    """
    if mut_str == 'WT' or pd.isna(mut_str):
        return wt_seq
    seq = list(wt_seq)
    # Try multiple separators
    for sep in ['/', ':', ',']:
        if sep in mut_str:
            parts = mut_str.split(sep)
            break
    else:
        parts = [mut_str]
    for p in parts:
        p = p.strip()
        r = parse_mut(p)
        if r is None:
            return None
        pos, new_aa = r
        if pos < 0 or pos >= len(seq):
            return None
        seq[pos] = new_aa
    return ''.join(seq)

sequences = []
errors = []
for i, row in df.iterrows():
    wt = WT.get(row['GFP type'])
    if wt is None:
        sequences.append(None)
        errors.append(f"unknown_type:{row['GFP type']}")
        continue
    seq = apply_mutations(wt, row['aaMutations'])
    if seq is None:
        sequences.append(None)
        errors.append(f"parse_fail:{row['aaMutations']}")
    else:
        sequences.append(seq)

print(f"  Done in {time.time()-t0:.1f}s")
n_err = sum(1 for s in sequences if s is None)
print(f"  Successful: {len(sequences) - n_err:,} / {len(sequences):,}")
print(f"  Errors: {n_err:,}")

if errors:
    from collections import Counter
    print("  Error breakdown:")
    for e, c in Counter(errors).most_common(10):
        print(f"    {c:>8,} {e}")

# ---- Build output dataframe ----
df_out = pd.DataFrame({
    'type': df['GFP type'].values,
    'mutations': df['aaMutations'].values,
    'brightness': df['Brightness'].values,
    'seq': sequences,
})
# Drop rows with failed seq
df_out = df_out.dropna(subset=['seq']).reset_index(drop=True)
print(f"\nFinal: {len(df_out):,} variants with valid sequences")

# ---- Save ----
import pickle
df_out.to_pickle(ROUND2 / "all_variants.pkl")
print(f"Saved to {ROUND2 / 'all_variants.pkl'}")

# Also save a quick CSV with just type + brightness + seq
df_small = df_out[['type', 'brightness', 'seq']].copy()
df_small.to_csv(ROUND2 / "all_variants.csv", index=False)
print(f"CSV (smaller): {ROUND2 / 'all_variants.csv'}")

# Quick stats per type
print("\nPer-type stats:")
print(df_out.groupby('type').agg(
    n=('brightness', 'count'),
    mean=('brightness', 'mean'),
    median=('brightness', 'median'),
    std=('brightness', 'std'),
    min=('brightness', 'min'),
    max=('brightness', 'max'),
).round(2))