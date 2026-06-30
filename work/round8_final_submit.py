"""
Round 8 Final: 合并所有候选 + 多样性 Top-6 + 提交
"""
import json, numpy as np, pandas as pd
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
R6 = ROOT / "work" / "round6"
R7 = ROOT / "work" / "round7"
R8 = ROOT / "work" / "round8"

# 加载所有候选
all_candidates = []

# 1. Round 6 (已有评分)
with open(R6 / "round6_full_ranking.json", encoding="utf-8") as f:
    all_candidates.extend(json.load(f))

# 2. Round 7 LMPNN
r7_path = R7 / "esmfold_lmpnn_r7.json"
if r7_path.exists():
    with open(r7_path, encoding="utf-8") as f:
        all_candidates.extend(json.load(f))

# 3. Round 7 ProteinMPNN
r7_mpnn = R7 / "esmfold_round7.json"
if r7_mpnn.exists():
    with open(r7_mpnn, encoding="utf-8") as f:
        all_candidates.extend(json.load(f))

# 4. Round 8 LMPNN (new)
r8_path = R8 / "esmfold_lmpnn_r8.json"
if r8_path.exists():
    with open(r8_path, encoding="utf-8") as f:
        all_candidates.extend(json.load(f))

print(f"Total: {len(all_candidates)} candidates")

# 新规评分
def compute(c):
    ptm = c.get("ptm", 0) or 0
    plddt = c.get("plddt_mean", 0) or 0
    chromo = c.get("plddt_chromo_region", plddt) or plddt
    nmuts = c.get("n_muts", 0) or 0
    return round(ptm * 0.50 + plddt / 100 * 0.30 + chromo / 100 * 0.20 - nmuts / 200 * 0.05, 4)

for c in all_candidates:
    c["sort_score"] = compute(c)
    c["pass_hard"] = (c.get("ptm") or 0) > 0.50 and (c.get("plddt_mean") or 0) > 50.0
    c["top2_ok"] = (c.get("ptm") or 0) > 0.70 and (c.get("plddt_mean") or 0) > 60

all_candidates.sort(key=lambda x: -x["sort_score"])

print(f"\nTop-20 全局排名:")
print(f"{'#':<4} {'Name':<35} {'Scaffold':<14} {'pLDDT':>6} {'pTM':>6} {'Muts':>4} {'Score':>8} {'Top2':>5}")
for i, c in enumerate(all_candidates[:20], 1):
    t = "✅" if c.get("top2_ok") else ""
    print(f"{i:<4} {c['name'][:35]:<35} {c.get('scaffold','?')[:14]:<14} "
          f"{c.get('plddt_mean',0):>6.1f} {c.get('ptm',0):>6.4f} "
          f"{c.get('n_muts',0):>4} {c['sort_score']:>8.4f} {t:>5}")

# ============================================================
# 多样性 Top-6 选择
# ============================================================
def cat(scaffold):
    s = str(scaffold)
    if "LMPNN" in s: return "LMPNN"
    if "MPNN" in s: return "MPNN"
    return "manual"

# 1. Top-2 达标候选池
top2_pool = [c for c in all_candidates if c.get("top2_ok")]
selected = top2_pool[:2]

# 2. 多样性补充
seen = set(c["seq"] for c in selected)
pool = [c for c in all_candidates if c["seq"] not in seen and (c.get("ptm") or 0) > 0.55]

for cat_name in ["MPNN", "LMPNN", "manual"]:
    for c in pool:
        if cat(c.get("scaffold","")) == cat_name:
            selected.append(c)
            seen.add(c["seq"])
            break

# 补到6条
if len(selected) < 6:
    for c in pool:
        if c["seq"] not in seen:
            selected.append(c)
            seen.add(c["seq"])
        if len(selected) >= 6:
            break

selected.sort(key=lambda x: -x["sort_score"])

print(f"\n{'='*120}")
print(f"FINAL Top-6 提交 (新规则)")
print(f"{'='*120}")
print(f"{'Seq':<5} {'Name':<35} {'Scaffold':<14} {'pLDDT':>6} {'Chromo':>6} "
      f"{'pTM':>6} {'Muts':>4} {'Score':>8} {'Status':<15}")
print("-" * 120)

for i, c in enumerate(selected, 1):
    ptm = c.get("ptm", 0) or 0
    plddt = c.get("plddt_mean", 0) or 0
    chromo = c.get("plddt_chromo_region", plddt) or plddt
    
    if ptm > 0.70 and plddt > 60:
        status = "✅ Top2达标"
    elif ptm > 0.55:
        status = "🟡 多样性"
    else:
        status = "🔴 不达标"
    
    print(f"{i:<5} {c['name'][:35]:<35} {c.get('scaffold','?')[:14]:<14} "
          f"{plddt:>6.1f} {chromo:>6.1f} {ptm:>6.4f} "
          f"{c.get('n_muts',0):>4} {c['sort_score']:>8.4f} {status:<15}")

# 提交文件
team_name = "YourTeamName"
sub = pd.DataFrame({
    "Team_Name": [team_name] * 6,
    "Seq_ID": list(range(1, 7)),
    "Sequence": [c["seq"] for c in selected]
})

# 合规检查
excl = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())
sub_path = R8 / "submission_round8_final.csv"
sub.to_csv(sub_path, index=False)
print(f"\n✅ 提交: {sub_path}")

# Top-2 达标数
n_top2 = sum(1 for c in selected if (c.get("ptm") or 0) > 0.70 and (c.get("plddt_mean") or 0) > 60)
print(f"Top2达标: {n_top2}/2")

# 保存 JSON
final = []
for i, c in enumerate(selected, 1):
    final.append({
        "Seq_ID": i, "name": c["name"],
        "scaffold": c.get("scaffold"), "n_muts": c.get("n_muts"),
        "plddt_mean": c.get("plddt_mean"),
        "plddt_chromo_region": c.get("plddt_chromo_region"),
        "ptm": c.get("ptm"), "sort_score": c["sort_score"],
        "seq": c["seq"],
    })
with open(R8 / "final_6_round8.json", "w", encoding="utf-8") as f:
    json.dump(final, f, indent=2, ensure_ascii=False)
print(f"✅ JSON: {R8 / 'final_6_round8.json'}")
