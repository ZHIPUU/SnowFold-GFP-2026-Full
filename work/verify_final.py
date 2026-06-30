"""Verify submission meets all competition requirements"""
import pandas as pd, json
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
PRE = ROOT / "预选序列"

sub = pd.read_csv(PRE / "submission_top6.csv")
excl = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())

print("=== 提交验证 ===")
print(f"序列数: {len(sub)}")
errors = []
for i, row in sub.iterrows():
    seq = row["Sequence"]
    if not seq.startswith("M"):
        errors.append(f"Seq {i+1}: 不以M开头")
    if not (220 <= len(seq) <= 250):
        errors.append(f"Seq {i+1}: 长度 {len(seq)} (需220-250)")
    bad = set(seq) - set("ACDEFGHIKLMNPQRSTVWY")
    if bad:
        errors.append(f"Seq {i+1}: 非法字符 {bad}")
    if seq in excl:
        errors.append(f"Seq {i+1}: 在排除列表中")

if errors:
    print("❌ 错误:")
    for e in errors:
        print(f"  {e}")
else:
    print("✅ 全部合规!")

with open(PRE / "final_top6.json") as f:
    final = json.load(f)
print(f"\nTop-6 JSON: {len(final)} 条")
for c in final:
    t2 = "✅" if c["ptm"] > 0.70 and c["plddt_mean"] > 60 else ""
    print(f"  Seq {c['Seq_ID']}: {c['name']} (pTM={c['ptm']}, pLDDT={c['plddt_mean']}) {t2}")
