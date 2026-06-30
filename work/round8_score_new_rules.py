"""
按新规则重新评分所有候选
===========================
新公式: 排序分 = pTM×0.50 + (pLDDT/100)×0.30 + (chromo/100)×0.20 - (muts/200)×0.05

新阈值:
  - 硬性通过: pTM>0.50, pLDDT>50.0, chromo>45.0
  - 正信号(候选池): pTM>0.65 且 chromo>55
  - 提交要求: Top 2 -> pTM>0.70 + pLDDT>60
  - 提交要求: 其余4 -> pTM>0.55

提交策略: 6条中只要求Top2达标, 其余4条展示多样性
"""
import json, numpy as np
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
R5 = ROOT / "work" / "round5"
R6 = ROOT / "work" / "round6"
R7 = ROOT / "work" / "round7"

# 加载所有已有候选
all_candidates = []

# Round 6 修正版已评分的
with open(R6 / "round6_full_ranking.json", encoding="utf-8") as f:
    all_candidates.extend(json.load(f))

# Round 7 LMPNN 新候选
lmpnn_path = R7 / "esmfold_lmpnn_r7.json"
if lmpnn_path.exists():
    with open(lmpnn_path, encoding="utf-8") as f:
        all_candidates.extend(json.load(f))

# Round 7 MPNN 新候选
mpnn_path = R7 / "esmfold_round7.json"
if mpnn_path.exists():
    with open(mpnn_path, encoding="utf-8") as f:
        all_candidates.extend(json.load(f))

print(f"总候选: {len(all_candidates)} 条")

# 新公式评分
def compute_sort_score(c):
    ptm = c.get("ptm", 0) or 0
    plddt = c.get("plddt_mean", 0) or 0
    chromo = c.get("plddt_chromo_region", plddt) or plddt
    n_muts = c.get("n_muts", 0) or 0
    
    score = (ptm * 0.50) + (plddt / 100 * 0.30) + (chromo / 100 * 0.20) - (n_muts / 200 * 0.05)
    return round(score, 4)

for c in all_candidates:
    c["sort_score"] = compute_sort_score(c)
    c["pass_hard"] = (c.get("ptm",0) or 0) > 0.50 and (c.get("plddt_mean",0) or 0) > 50.0 and (c.get("plddt_chromo_region", c.get("plddt_mean",0)) or 0) > 45.0
    c["pass_positive"] = (c.get("ptm",0) or 0) > 0.65 and (c.get("plddt_chromo_region", c.get("plddt_mean",0)) or 0) > 55

# 排序
all_candidates.sort(key=lambda x: -x["sort_score"])

# 检查提交达标数
top2_ok = [c for c in all_candidates if (c.get("ptm",0) or 0) > 0.70 and (c.get("plddt_mean",0) or 0) > 60]
remaining_ok = [c for c in all_candidates if (c.get("ptm",0) or 0) > 0.55 and not ((c.get("ptm",0) or 0) > 0.70 and (c.get("plddt_mean",0) or 0) > 60)]

print(f"\nTop 2 达标 (pTM>0.70 + pLDDT>60): {len(top2_ok)} 条")
print(f"剩余达标 (pTM>0.55): {len(remaining_ok)} 条")
print(f"硬性通过 (pTM>0.50 + pLDDT>50 + chromo>45): {sum(1 for c in all_candidates if c['pass_hard'])} 条")
print(f"正信号 (pTM>0.65 + chromo>55): {sum(1 for c in all_candidates if c['pass_positive'])} 条")

print(f"\n{'='*130}")
print(f"新规排序 Top-30")
print(f"{'='*130}")
print(f"{'#':<4} {'Name':<35} {'Scaffold':<14} {'pLDDT':>6} {'Chromo':>6} {'pTM':>6} "
      f"{'Muts':>4} {'SortScore':>9} {'Hard':>5} {'Pos':>5}")
print("-" * 130)
for i, c in enumerate(all_candidates[:30], 1):
    ptm = c.get("ptm", 0) or 0
    plddt = c.get("plddt_mean", 0) or 0
    chromo = c.get("plddt_chromo_region", plddt) or plddt
    hard = "✅" if c["pass_hard"] else "❌"
    pos = "✅" if c["pass_positive"] else "❌"
    print(f"{i:<4} {c['name'][:35]:<35} {c.get('scaffold','?')[:14]:<14} "
          f"{plddt:>6.1f} {chromo:>6.1f} {ptm:>6.4f} "
          f"{c.get('n_muts',0):>4} {c['sort_score']:>9.4f} {hard:>5} {pos:>5}")

# ============================================================
# 多样性 Top-6 选择 (新规则)
# ============================================================
def cat(scaffold):
    s = str(scaffold)
    if "MPNN" in s and "LMPNN" not in s: return "MPNN"
    elif "LMPNN" in s: return "LMPNN"
    elif s in ("sfGFP", "avGFP", "amacGFP", "cgreGFP"): return "manual"
    return s

# 1. 先确保 Top 2 达标 (pTM>0.70 + pLDDT>60)
top2_pool = [c for c in all_candidates if (c.get("ptm",0) or 0) > 0.70 and (c.get("plddt_mean",0) or 0) > 60]
selected = top2_pool[:2]  # 最多取2条

# 2. 补足到6条: 多样性 + pTM>0.55
seen_seqs = set(c["seq"] for c in selected)
pool = [c for c in all_candidates if (c.get("ptm",0) or 0) > 0.55 and c["seq"] not in seen_seqs]

# 按类别多样性
for cat_name in ["MPNN", "LMPNN", "manual"]:
    for c in pool:
        if c["seq"] not in seen_seqs and cat(c.get("scaffold","")) == cat_name:
            selected.append(c)
            seen_seqs.add(c["seq"])
        if len(selected) >= 6:
            break
    if len(selected) >= 6:
        break

# 如果还不够, 从 pool 补
if len(selected) < 6:
    for c in pool:
        if c["seq"] not in seen_seqs:
            selected.append(c)
            seen_seqs.add(c["seq"])
        if len(selected) >= 6:
            break

selected.sort(key=lambda x: -x["sort_score"])

print(f"\n{'='*120}")
print(f"🎯 新规 Top-6 提交")
print(f"{'='*120}")
print(f"{'Seq':<5} {'Name':<35} {'Scaffold':<14} {'pLDDT':>6} {'Chromo':>6} "
      f"{'pTM':>6} {'Muts':>4} {'Score':>9} {'Qualify':<12}")
print("-" * 120)

for i, c in enumerate(selected, 1):
    ptm = c.get("ptm", 0) or 0
    plddt = c.get("plddt_mean", 0) or 0
    chromo = c.get("plddt_chromo_region", plddt) or plddt
    
    if ptm > 0.70 and plddt > 60:
        qualify = "✅ Top2达标"
    elif ptm > 0.55:
        qualify = "🟡 多样性"
    else:
        qualify = "🔴 不达标"
    
    print(f"{i:<5} {c['name'][:35]:<35} {c.get('scaffold','?')[:14]:<14} "
          f"{plddt:>6.1f} {chromo:>6.1f} {ptm:>6.4f} "
          f"{c.get('n_muts',0):>4} {c['sort_score']:>9.4f} {qualify:<12}")

# 提交文件
import pandas as pd
sub = pd.DataFrame({
    "Team_Name": ["YourTeamName"] * 6,
    "Seq_ID": list(range(1, 7)),
    "Sequence": [c["seq"] for c in selected]
})

# 合规
excl = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())
for s in sub["Sequence"]:
    assert s.startswith("M")
    assert 220 <= len(s) <= 250
    assert not (set(s) - set("ACDEFGHIKLMNPQRSTVWY"))
    assert s not in excl, "Excluded!"

sub_path = R6 / "submission_new_rules.csv"
sub.to_csv(sub_path, index=False)
print(f"\n提交: {sub_path} | 合规: ✅")

# 保存 JSON
final = []
for i, c in enumerate(selected, 1):
    final.append({
        "Seq_ID": i, "name": c["name"], "scaffold": c.get("scaffold"),
        "n_muts": c.get("n_muts"), "plddt_mean": c.get("plddt_mean"),
        "plddt_chromo_region": c.get("plddt_chromo_region"),
        "ptm": c.get("ptm"), "sort_score": c["sort_score"],
        "seq": c["seq"],
    })
with open(R6 / "final_6_new_rules.json", "w", encoding="utf-8") as f:
    json.dump(final, f, indent=2, ensure_ascii=False)

print(f"\nTop-2 达标: {sum(1 for c in selected if (c.get('ptm',0) or 0) > 0.70 and (c.get('plddt_mean',0) or 0) > 60)} / 2")
