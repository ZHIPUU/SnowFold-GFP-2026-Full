"""
Round 3 最终提交生成
从 8 条候选中选择最优 6 条，基于 pLDDT + 多样性
"""
import json, pandas as pd
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")

# 加载 ESMFold 结果
with open(ROOT / "work/round3/esmfold_results.json") as f:
    fold_results = json.load(f)

# 加载候选
with open(ROOT / "work/round3/candidates_round3.json") as f:
    candidates = json.load(f)

# 合并
cand_map = {c["name"]: c for c in candidates}
for r in fold_results:
    r["notes"] = cand_map[r["name"]]["notes"]
    r["seq"] = cand_map[r["name"]]["seq"]

# 过滤 pLDDT >= 70
good = [r for r in fold_results if r["plddt_mean"] >= 70]
print(f"pLDDT >= 70 的候选: {len(good)}/{len(fold_results)}")
for r in fold_results:
    status = "✓" if r["plddt_mean"] >= 70 else "✗"
    print(f"  {status} {r['plddt_mean']:5.1f} {r['name']}")

# 分类 scaffold
from collections import defaultdict
by_scaffold = defaultdict(list)
for r in good:
    for s in ["avGFP", "sfGFP", "amacGFP", "cgreGFP", "mBaoJin"]:
        if r["name"].startswith(s):
            by_scaffold[s].append(r)
            break

print(f"\n各 scaffold 有效候选:")
for s, items in by_scaffold.items():
    print(f"  {s}: {len(items)}")

# 选 6 条策略:
# - 每 scaffold 至少 1 条
# - pLDDT 降序选，优先高 pLDDT
# - 如果某 scaffold 无合格候选，用其他补

selected = []
scaffolds_selected = set()

# 第一轮: 每 scaffold 选最佳
for s, items in by_scaffold.items():
    if items:
        best = max(items, key=lambda x: x["plddt_mean"])
        selected.append(best)
        scaffolds_selected.add(s)

# 第二轮: 如果不足 6 条，按 pLDDT 从剩余中补
remaining = [r for r in good if r not in selected]
remaining.sort(key=lambda x: x["plddt_mean"], reverse=True)
while len(selected) < 6 and remaining:
    selected.append(remaining.pop(0))

# 如果还不够，降到 plddt >= 65
if len(selected) < 6:
    fallback = [r for r in fold_results if r["plddt_mean"] >= 65 and r not in selected]
    fallback.sort(key=lambda x: x["plddt_mean"], reverse=True)
    selected.extend(fallback[:6 - len(selected)])

print(f"\n最终选择 {len(selected)} 条:")
for i, r in enumerate(selected, 1):
    print(f"  Seq {i}: {r['name']} (pLDDT={r['plddt_mean']:.1f})")
    print(f"    {r['notes']}")

# 生成提交 CSV
sub = pd.DataFrame({
    "Team_Name": ["YourTeamName"] * len(selected),
    "Seq_ID": list(range(1, len(selected) + 1)),
    "Sequence": [r["seq"] for r in selected]
})

out_path = ROOT / "submission_yourteamname.csv"
sub.to_csv(out_path, index=False)
print(f"\n提交已保存到 {out_path}")
print(f"共 {len(selected)} 条序列")

# 快速验证
print("\n=== 快速验证 ===")
for _, row in sub.iterrows():
    s = row["Sequence"]
    chromo = "?"
    for t in ["TYG", "SYG", "GYG"]:
        if t in s: chromo = t; break
    ok = s.startswith("M") and 220 <= len(s) <= 250
    print(f"  Seq {row['Seq_ID']}: len={len(s)} M-start={s[0]=='M'} chromo={chromo} {'✓' if ok else '✗'}")
