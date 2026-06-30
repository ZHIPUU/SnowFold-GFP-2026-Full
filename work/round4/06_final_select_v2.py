"""
Round 4 优化版 Top-6 最终选择
==================================
基于多样性诊断 + ESMFold v2 评估
策略 (按 Best Top-1 + 多样性):
  1 条 X4_avGFP_sfGFP_I152S_Q69L (新王者, pLDDT最高, 12mut)
  1 条 C1_sfGFP_I152S (Round 3 验证保险)
  1 条 mBaoJin 系列 (热稳爆款, Tm 92°C)
  1 条 amacGFP 系列 (跨骨架多样)
  1 条 avGFP-only 不同变体 (再多样)
  1 条 B3 综合 (S30R + sf:acid)
"""
import json, pandas as pd
from pathlib import Path
from collections import defaultdict

ROOT = Path(r"D:\生信\2026Protein Design")
OUT = ROOT / "work" / "round4"

with open(OUT / "esmfold_round4_v2.json", encoding="utf-8") as f:
    results = json.load(f)

# 删除明显失败的 cgreGFP (pLDDT ~30, 结构有问题)
results = [r for r in results if r["scaffold"] != "cgreGFP"]
print(f"剔除 cgreGFP 后: {len(results)} 条")

def score(r):
    plddt = r["plddt_mean"]
    chromo = r["plddt_chromo_region"]
    ptm = r["ptm"] or 0.0
    n_mut = r["n_muts"]
    tm = r["expected_tm"]

    plddt_score = max(0, min(1, (plddt - 35) / 20))
    chromo_score = max(0, min(1, (chromo - 30) / 35))
    ptm_score = max(0, min(1, (ptm - 0.3) / 0.3))
    tm_score = max(0, min(1, (tm - 70) / 25))
    if n_mut <= 5:
        brightness_score = 1.0
    elif n_mut <= 10:
        brightness_score = 1.0 - (n_mut - 5) * 0.05
    else:
        brightness_score = max(0, 0.75 - (n_mut - 10) * 0.15)

    weights = {"plddt": 1.5, "chromo": 2.0, "ptm": 1.0, "tm": 2.5, "brightness": 2.0}
    total = (weights["plddt"] * plddt_score + weights["chromo"] * chromo_score +
             weights["ptm"] * ptm_score + weights["tm"] * tm_score +
             weights["brightness"] * brightness_score)
    return round(total / sum(weights.values()) * 10, 2)

for r in results:
    r["score"] = score(r)

sorted_r = sorted(results, key=lambda x: -x["score"])

print("\n" + "=" * 110)
print("Top-15 候选 (按综合得分):")
print("=" * 110)
print(f"{'#':<3} {'name':<32} {'scaffold':<10} {'mut':>3} {'pLDDT':>6} {'cb':>6} {'pTM':>5} {'Tm':>3} {'score':>5}")
print("-" * 110)
for i, r in enumerate(sorted_r[:15], 1):
    print(f"{i:<3} {r['name']:<32} {r['scaffold']:<10} {r['n_muts']:>3} "
          f"{r['plddt_mean']:>6.2f} {r['plddt_chromo_region']:>6.2f} {r['ptm']:>5.3f} "
          f"{r['expected_tm']:>3} {r['score']:>5.2f}")

# ============================================================
# 多样性优化选择 (Diverse Top-6)
# 策略: 每个骨架最多 2 条, 强制至少 3 种骨架
# ============================================================
print("\n" + "=" * 70)
print("多样性优化 Top-6 (每个骨架最多 2 条):")
print("=" * 70)

selected = []
seen_seqs = set()
scaffold_count = defaultdict(int)

# Phase 1: 每个骨架取最佳 1 条 (强制多样性)
print("\nPhase 1 - 每个骨架取最佳 (强制多样性):")
for scaffold in ["sfGFP", "avGFP", "amacGFP", "mBaoJin"]:
    for r in sorted_r:
        if r["scaffold"] == scaffold and r["seq"] not in seen_seqs:
            selected.append(r)
            seen_seqs.add(r["seq"])
            scaffold_count[scaffold] += 1
            print(f"  ✓ {r['name']:<32} ({scaffold}) score={r['score']}")
            break

# Phase 2: 补足到 6, 优先高分 + 不同序列, 每骨架≤2
print(f"\nPhase 2 - 补足到 6 (每骨架最多 2):")
for r in sorted_r:
    if len(selected) >= 6: break
    if r["seq"] in seen_seqs: continue
    if scaffold_count[r["scaffold"]] >= 2: continue  # 每骨架最多 2 条
    selected.append(r)
    seen_seqs.add(r["seq"])
    scaffold_count[r["scaffold"]] += 1
    print(f"  ➕ {r['name']:<32} ({r['scaffold']}) score={r['score']}")

# 按分数重排 (Seq_ID 1 = 最高分)
final_6 = sorted(selected[:6], key=lambda x: -x["score"])

print("\n" + "=" * 110)
print("🎉 Round 4 优化版 最终 6 条:")
print("=" * 110)
print(f"{'Seq_ID':<6} {'name':<32} {'scaffold':<10} {'mut':>3} {'pLDDT':>6} {'cb':>6} {'Tm':>3} {'score':>5}")
print("-" * 110)
for i, r in enumerate(final_6, 1):
    print(f"{i:<6} {r['name']:<32} {r['scaffold']:<10} {r['n_muts']:>3} "
          f"{r['plddt_mean']:>6.2f} {r['plddt_chromo_region']:>6.2f} {r['expected_tm']:>3} {r['score']:>5.2f}")

# 多样性检查
print("\n多样性诊断:")
sf_count = defaultdict(int)
for r in final_6:
    sf_count[r["scaffold"]] += 1
for sc, n in sf_count.items():
    print(f"  {sc}: {n} 条")

# 汉明距离矩阵
print("\n汉明距离矩阵 (越大越好, 多样性):")
print(f"{'':<6}", end="")
for i in range(6): print(f"{i+1:>5}", end="")
print()
for i, c1 in enumerate(final_6):
    print(f"Seq{i+1:<3}", end="  ")
    for j, c2 in enumerate(final_6):
        if i == j:
            print(f"{'-':>5}", end="")
        else:
            min_len = min(len(c1["seq"]), len(c2["seq"]))
            d = sum(1 for a, b in zip(c1["seq"][:min_len], c2["seq"][:min_len]) if a != b)
            d += abs(len(c1["seq"]) - len(c2["seq"]))
            print(f"{d:>5}", end="")
    print()

# 保存
output_data = []
for i, r in enumerate(final_6, 1):
    output_data.append({
        "Seq_ID": i, "name": r["name"], "role": r["role"],
        "scaffold": r["scaffold"], "n_muts": r["n_muts"],
        "expected_tm": r["expected_tm"],
        "plddt_mean": r["plddt_mean"], "plddt_chromo_region": r["plddt_chromo_region"],
        "ptm": r["ptm"], "round4_score": r["score"],
        "notes": r["notes"], "seq": r["seq"],
    })

with open(OUT / "final_6_round4_v2.json", "w", encoding="utf-8") as f:
    json.dump(output_data, f, indent=2, ensure_ascii=False)

sub = pd.DataFrame({
    "Team_Name": ["YourTeamName"] * 6,
    "Seq_ID": list(range(1, 7)),
    "Sequence": [r["seq"] for r in final_6]
})
sub_path = OUT / "submission_round4_v2.csv"
sub.to_csv(sub_path, index=False)

# 合规性检查
print("\n" + "=" * 70)
print("提交合规性检查:")
print("=" * 70)
excl_seqs = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())
all_ok = True
for _, row in sub.iterrows():
    s = row["Sequence"]
    issues = []
    if not s.startswith("M"): issues.append("not-M")
    if not (220 <= len(s) <= 250): issues.append(f"len={len(s)}")
    if set(s) - set("ACDEFGHIKLMNPQRSTVWY"): issues.append("bad-AA")
    if s in excl_seqs: issues.append("EXCLUDED!")
    status = "✅ OK" if not issues else "❌ " + ", ".join(issues)
    print(f"  Seq {row['Seq_ID']}: len={len(s)}  {status}")
    if issues: all_ok = False

print(f"\n{'🎉 全部通过!' if all_ok else '⚠️ 有问题'}")
print(f"\n📁 提交文件: {sub_path}")
