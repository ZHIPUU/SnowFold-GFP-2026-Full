"""
Phase 1.6: 筛选 Top-50 跨母体保守候选
================================================
策略:
  - 只保留 2-4 突变 (加性模型在低突变数下更可靠)
  - 跨母体均匀分配 (~12 per type)
  - 去重序列
  - 输出供 Phase 2 物理验证
"""
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

def parse_mut(s):
    return []
sys.modules['__main__'].parse_mut = parse_mut

WORK = Path(r"D:\生信\2026Protein Design\work")
PHASE1 = WORK / "phase1"

# 读取全候选
df = pd.read_csv(PHASE1 / "top_candidates.csv")
print(f"Total candidates: {len(df)}")
print(f"Distribution:")
print(df.groupby(["type", "n_mut"]).size().unstack(fill_value=0))

# 只保留 2-4 突变
df_filter = df[df["n_mut"].isin([2, 3, 4])].copy()
print(f"\nAfter filter 2-4 mut: {len(df_filter)}")

# 按 type + n_mut 分组, 取 top-N
TOP_PER_GROUP = 12
selected = []
for t in ["avGFP", "amacGFP", "cgreGFP", "ppluGFP"]:
    for n in [2, 3, 4]:
        sub = df_filter[(df_filter["type"] == t) & (df_filter["n_mut"] == n)]
        sub = sub.sort_values("pred_brightness", ascending=False).head(TOP_PER_GROUP)
        selected.append(sub)

final = pd.concat(selected, ignore_index=True)
print(f"\nFinal selected: {len(final)}")

# 去重 seq
final = final.drop_duplicates(subset="seq").reset_index(drop=True)
print(f"After dedupe: {len(final)}")

# 输出
final.to_csv(PHASE1 / "top50_candidates.csv", index=False)
print(f"Saved to {PHASE1 / 'top50_candidates.csv'}")

# 输出分布
print("\n=== Final candidates by type/n_mut ===")
print(final.groupby(["type", "n_mut"]).size().unstack(fill_value=0))

print("\n=== Final top 20 ===")
for _, r in final.head(20).iterrows():
    print(f"  {r['type']:8s} {r['n_mut']}-mut: pred={r['pred_brightness']:.3f} (rel={r['pred_relative']:.2f}×)  {r['mut_str']}")