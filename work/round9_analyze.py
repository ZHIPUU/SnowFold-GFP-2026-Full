"""Analyze Round 9 results and produce final submission"""
import json
from pathlib import Path

R9 = Path(r"D:\生信\2026Protein Design\work\round9")
R8 = Path(r"D:\生信\2026Protein Design\work\round8")

# Load
with open(R9 / "esmfold_r9.json") as f:
    data = json.load(f)

# Best <30 muts
best = sorted([v for v in data if v["n_muts"] < 30], key=lambda x: -x["ptm"])
print("Best low-mutation variants (<30 muts):")
for v in best[:5]:
    print(f"  {v['name']}: muts={v['n_muts']}, pLDDT={v['plddt_mean']:.1f}, pTM={v['ptm']:.4f}")

# Score comparison with current Seq 6
b = best[0]
new_score = round(b["ptm"] * 0.50 + b["plddt_mean"]/100 * 0.30 + b.get("plddt_chromo", b["plddt_mean"])/100 * 0.20 - b["n_muts"]/200 * 0.05, 4)

with open(R8 / "final_6_round8.json") as f:
    final = json.load(f)
seq6 = final[-1]
old_score = seq6["sort_score"]

print(f"\n替换建议:")
print(f"  当前 Seq 6: {seq6['name']}, pTM={seq6['ptm']}, score={old_score}")
print(f"  候选替换:   {b['name']}, muts={b['n_muts']}, pTM={b['ptm']:.4f}, score={new_score}")
print(f"  提升:       score {old_score} -> {new_score} ({'+' if new_score > old_score else ''}{new_score - old_score:.4f})")
