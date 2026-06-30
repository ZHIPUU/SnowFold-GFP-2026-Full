"""
Round 4 v3 终极版 Top-6 选择
==================================
合并所有评估结果:
  - 41 条手工设计 (sfGFP / avGFP / amacGFP / mBaoJin / cgreGFP)
  - 12 条 ProteinMPNN de novo (sfGFP_MPNN)
合计: 53 条候选

策略 (Best Top-1 + 多样性):
  - 1 条 MPNN de novo 王者 (pLDDT 68 最高)
  - 1 条 sfGFP 保险 (Round 3 验证)
  - 1 条 avGFP+sfGFP10+多机制叠加 (X4 王者)
  - 1 条 amacGFP / 跨骨架
  - 1 条 mBaoJin 高Tm爆款
  - 1 条 htFuncLib 综合
"""
import json, pandas as pd
from pathlib import Path
from collections import defaultdict

ROOT = Path(r"D:\生信\2026Protein Design")
OUT = ROOT / "work" / "round4"

# 加载所有评估结果
with open(OUT / "esmfold_round4_v3.json", encoding="utf-8") as f:
    handcraft = json.load(f)
with open(OUT / "esmfold_mpnn.json", encoding="utf-8") as f:
    mpnn = json.load(f)

# 给 MPNN 候选补 n_muts (相对 sfGFP) 已有
# 给 handcraft 加 mpnn_recovery 字段防止报错
for r in handcraft:
    r.setdefault("mpnn_recovery", None)
    r.setdefault("mpnn_score", None)

all_results = handcraft + mpnn
print(f"总评估候选: {len(all_results)}")

# 剔除明显失败的 cgreGFP (pLDDT < 35)
all_results = [r for r in all_results if r["plddt_mean"] >= 35]
print(f"过滤后 (pLDDT>=35): {len(all_results)}")

def score(r):
    plddt = r["plddt_mean"]
    chromo = r["plddt_chromo_region"]
    ptm = r["ptm"] or 0.0
    n_mut = r["n_muts"]
    tm = r["expected_tm"]
    plddt_s = max(0, min(1, (plddt - 35) / 30))  # 35-65 区间
    chromo_s = max(0, min(1, (chromo - 30) / 40))  # 30-70
    ptm_s = max(0, min(1, (ptm - 0.3) / 0.5))  # 0.3-0.8
    tm_s = max(0, min(1, (tm - 70) / 25))
    if n_mut <= 5: br = 1.0
    elif n_mut <= 12: br = 1.0 - (n_mut - 5) * 0.03
    else: br = max(0, 0.79 - (n_mut - 12) * 0.05)
    weights = {"plddt": 2.0, "chromo": 2.0, "ptm": 1.5, "tm": 2.5, "brightness": 2.0}
    total = (weights["plddt"]*plddt_s + weights["chromo"]*chromo_s + weights["ptm"]*ptm_s 
             + weights["tm"]*tm_s + weights["brightness"]*br)
    return round(total / sum(weights.values()) * 10, 2)

for r in all_results:
    r["round4_score"] = score(r)

# 排序
sorted_r = sorted(all_results, key=lambda x: -x["round4_score"])
print("\n" + "=" * 115)
print("Top-20 综合排序:")
print("=" * 115)
print(f"{'#':<3} {'name':<28} {'scaffold':<12} {'mut':>3} {'pLDDT':>5} {'cb':>5} {'pTM':>5} {'Tm':>3} {'score':>5}")
print("-" * 115)
for i, r in enumerate(sorted_r[:20], 1):
    print(f"{i:<3} {r['name'][:28]:<28} {r['scaffold'][:12]:<12} {r['n_muts']:>3} "
          f"{r['plddt_mean']:>5.1f} {r['plddt_chromo_region']:>5.1f} {r['ptm'] or 0:>5.3f} "
          f"{r['expected_tm']:>3} {r['round4_score']:>5.2f}")

# ============================================================
# 多样性优化 Top-6
# ============================================================
print("\n" + "=" * 70)
print("Top-6 多样性优化 (每骨架≤2条):")
print("=" * 70)

scaffold_priority = ["sfGFP_MPNN", "avGFP", "sfGFP", "amacGFP", "mBaoJin"]

selected = []
seen_seqs = set()
scaffold_count = defaultdict(int)

# Phase 1: 每个骨架取最高分 1 条
for sc in scaffold_priority:
    for r in sorted_r:
        if r["scaffold"] == sc and r["seq"] not in seen_seqs:
            selected.append(r)
            seen_seqs.add(r["seq"])
            scaffold_count[sc] += 1
            print(f"  ✓ {r['name'][:28]:<28} ({sc}) score={r['round4_score']}")
            break

# Phase 2: 补足到 6 (每骨架≤2)
for r in sorted_r:
    if len(selected) >= 6: break
    if r["seq"] in seen_seqs: continue
    if scaffold_count[r["scaffold"]] >= 2: continue
    selected.append(r)
    seen_seqs.add(r["seq"])
    scaffold_count[r["scaffold"]] += 1
    print(f"  ➕ {r['name'][:28]:<28} ({r['scaffold']}) score={r['round4_score']}")

# 按 score 重排
final_6 = sorted(selected[:6], key=lambda x: -x["round4_score"])

print("\n" + "=" * 115)
print("🎉 Round 4 v3 终极版 Top-6:")
print("=" * 115)
print(f"{'Seq':<4} {'name':<28} {'scaffold':<12} {'mut':>3} {'pLDDT':>5} {'cb':>5} {'Tm':>3} {'score':>5}")
print("-" * 115)
for i, r in enumerate(final_6, 1):
    print(f"{i:<4} {r['name'][:28]:<28} {r['scaffold'][:12]:<12} {r['n_muts']:>3} "
          f"{r['plddt_mean']:>5.1f} {r['plddt_chromo_region']:>5.1f} "
          f"{r['expected_tm']:>3} {r['round4_score']:>5.2f}")

# 保存
output_data = []
for i, r in enumerate(final_6, 1):
    output_data.append({
        "Seq_ID": i, "name": r["name"], "scaffold": r["scaffold"],
        "n_muts": r["n_muts"], "expected_tm": r["expected_tm"],
        "plddt_mean": r["plddt_mean"], "plddt_chromo_region": r["plddt_chromo_region"],
        "ptm": r["ptm"], "round4_score": r["round4_score"],
        "notes": r.get("notes",""), "seq": r["seq"],
        "mpnn_recovery": r.get("mpnn_recovery"),
        "mpnn_score": r.get("mpnn_score"),
    })

with open(OUT / "final_6_round4_v3.json", "w", encoding="utf-8") as f:
    json.dump(output_data, f, indent=2, ensure_ascii=False)

sub = pd.DataFrame({
    "Team_Name": ["YourTeamName"] * 6,
    "Seq_ID": list(range(1, 7)),
    "Sequence": [r["seq"] for r in final_6]
})
sub_path = OUT / "submission_round4_v3.csv"
sub.to_csv(sub_path, index=False)

# 合规性
print("\n📋 合规性最终检查:")
excl_seqs = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())
all_ok = True
for _, row in sub.iterrows():
    s = row["Sequence"]
    issues = []
    if not s.startswith("M"): issues.append("not-M")
    if not (220 <= len(s) <= 250): issues.append(f"len={len(s)}")
    if set(s) - set("ACDEFGHIKLMNPQRSTVWY"): issues.append("bad-AA")
    if s in excl_seqs: issues.append("EXCL!")
    chromo = next((c for c in ["TYG","SYG","GYG"] if c in s), "?")
    status = "✅" if not issues else "❌"
    print(f"  Seq {row['Seq_ID']}: len={len(s)} chromo={chromo} {status}" + (" " + ",".join(issues) if issues else ""))
    if issues: all_ok = False

print(f"\n{'🎉 全部通过!' if all_ok else '⚠️ 有问题'}")
print(f"\n📁 提交文件: {sub_path}")
