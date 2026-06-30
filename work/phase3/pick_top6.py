"""从 phase2_scored.csv 选 top 6 跨母体候选"""
import sys
from pathlib import Path

import pandas as pd

WORK = Path(r"D:\生信\2026Protein Design\work")
PHASE2 = WORK / "phase2"
PHASE3 = WORK / "phase3"

cands = pd.read_csv(PHASE2 / "phase2_scored.csv")
print(f"Total scored: {len(cands)}")
print(cands.groupby("type").size())

# 看每个 type 的 best
for t in ["avGFP", "amacGFP", "cgreGFP", "ppluGFP"]:
    sub = cands[cands["type"] == t].sort_values("combined", ascending=False)
    print(f"\n=== {t} top 3 by combined ===")
    for _, r in sub.head(3).iterrows():
        print(f"  {r['n_mut']}-mut: pred={r['pred_brightness']:.3f}, esm_d={r['esm_delta_pll']:+.4f}, combined={r['combined']:.3f}")
        print(f"    mut: {r['mut_str']}")
        print(f"    seq_len: {len(r['seq'])}")