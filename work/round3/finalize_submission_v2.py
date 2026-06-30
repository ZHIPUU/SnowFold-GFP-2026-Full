"""
Round 3 最终提交生成 v2
按 pLDDT 相对排序选 6 条 + 多 scaffold 覆盖
"""
import json, pandas as pd
from pathlib import Path
from collections import defaultdict

ROOT = Path(r"D:\生信\2026Protein Design")

# 加载数据
with open(ROOT / "work/round3/esmfold_results.json") as f:
    fold_results = json.load(f)
with open(ROOT / "work/round3/candidates_round3.json") as f:
    candidates = json.load(f)

cand_map = {c["name"]: c for c in candidates}
for r in fold_results:
    r["notes"] = cand_map[r["name"]]["notes"]
    r["seq"] = cand_map[r["name"]]["seq"]

# 按 pLDDT 降序排序
fold_results.sort(key=lambda x: x["plddt_mean"], reverse=True)

print("=" * 70)
print("所有候选 pLDDT 排序")
print("=" * 70)
for r in fold_results:
    print(f"  pLDDT={r['plddt_mean']:5.1f}  {r['name']}")

# 分类 scaffold
def get_scaffold(name):
    for s in ["avGFP", "sfGFP", "amacGFP", "cgreGFP", "mBaoJin"]:
        if name.startswith(s):
            return s
    return "?"

by_scaffold = defaultdict(list)
for r in fold_results:
    by_scaffold[get_scaffold(r["name"])].append(r)

print("\nscaffold 分布:")
for s, items in sorted(by_scaffold.items(), key=lambda x: -len(x[1])):
    print(f"  {s}: {len(items)} 最佳={items[0]['plddt_mean']:.1f}")

# 选 6 条策略:
# 1. avGFP 系列取最佳 2-3 条 (sfGFP10 + sfGFP4core+S30R + sfGFP4core)
# 2. sfGFP + I152S (sfGFP 系列代表)
# 3. amacGFP+sfGFP5 (amacGFP 代表)
# 4. mBaoJin+D173N 或 cgreGFP (跨 scaffold)
# 排除 pLDDT 最低的 2 条 (cgreGFP 30.9, mBaoJin 38.8)

selected_names = [
    "sfGFP+I152S",           # 48.3 - 最佳, sfGFP 风格
    "amacGFP+sfGFP5",        # 47.5 - amacGFP 跨母体
    "avGFP+sfGFP10",         # 45.5 - avGFP 最强
    "avGFP+sfGFP4core+S30R", # 45.3 - 稳定增强
    "avGFP+sfGFP4core",      # 44.7 - 保守基础
    "avGFP+sfGFP4core+I152S", # 42.6 - 适度
]

selected = [r for r in fold_results if r["name"] in selected_names]
selected.sort(key=lambda x: [n for n in selected_names].index(x["name"]))

# 验证所有
print("\n" + "=" * 70)
print(f"最终选择 {len(selected)} 条")
print("=" * 70)
for i, r in enumerate(selected, 1):
    print(f"  Seq {i}: {r['name']} (pLDDT={r['plddt_mean']:.1f})")
    print(f"    {r['notes'][:90]}")

# 生成提交 CSV
sub = pd.DataFrame({
    "Team_Name": ["YourTeamName"] * len(selected),
    "Seq_ID": list(range(1, len(selected) + 1)),
    "Sequence": [r["seq"] for r in selected]
})

out_path = ROOT / "submission_yourteamname.csv"
sub.to_csv(out_path, index=False)
print(f"\n提交已保存到: {out_path}")
print(f"共 {len(selected)} 条序列")

# 快速验证
print("\n=== 快速验证 ===")
excl = pd.read_csv(ROOT / "Exclusion_List.csv")
excl_seqs = set(excl["Sequence"].astype(str).str.strip())

for _, row in sub.iterrows():
    s = row["Sequence"]
    chromo = "?"
    for t in ["TYG", "SYG", "GYG"]:
        if t in s: chromo = t; break
    in_excl = s in excl_seqs
    issues = []
    if not s.startswith("M"): issues.append("M?")
    if len(s) < 220 or len(s) > 250: issues.append(f"len={len(s)}")
    ok = "✓" if not issues else "✗"
    excl_ok = "✓" if not in_excl else "❌命中排除!"
    print(f"  Seq {row['Seq_ID']}: len={len(s)} chromo={chromo} 合规={ok} 排除={excl_ok}")

print(f"\n提交文件: {out_path}")
print(f"文件大小: {out_path.stat().st_size} bytes")