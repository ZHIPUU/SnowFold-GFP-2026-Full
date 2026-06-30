"""
Round 4 Step 3: 综合打分 + Top-6 选择
==================================
评分维度 (针对 Best Top-1 规则):
  1. pLDDT_score: ESMFold 整体置信度
  2. chromo_score: chromophore 区域 pLDDT (折叠正确性更直接指标)
  3. ptm_score: 全局拓扑置信度
  4. tm_score: 预期 Tm 文献先验 (热稳奖关键)
  5. brightness_score: 突变数惩罚 + sfGFP 兼容性
  6. role_diversity: 6 条要覆盖不同角色

策略 (Best Top-1):
  - 至少 1 条"高Tm爆款"(mBaoJin) - 综合分上限高
  - 至少 1 条"保险条"(sfGFP+I152S) - Round 3 验证
  - 2-3 条"htFuncLib风格"(综合分平衡)
  - 1 条"教科书"(avGFP+sfGFP10)
"""
import json, pandas as pd
from pathlib import Path
from collections import defaultdict

ROOT = Path(r"D:\生信\2026Protein Design")
OUT = ROOT / "work" / "round4"

with open(OUT / "esmfold_round4.json", encoding="utf-8") as f:
    results = json.load(f)

# ============================================================
# 评分函数
# ============================================================
def score_candidate(r):
    """
    综合得分 (0-10, 越高越好)
    专为 Best Top-1 优化: 上限 + 稳定性
    """
    plddt = r["plddt_mean"]
    chromo = r["plddt_chromo_region"]
    ptm = r["ptm"] or 0.0
    n_mut = r["n_muts"]
    tm = r["expected_tm"]

    # 1. pLDDT 折叠置信度 (规范化到 0-1)
    #    GFP 上 ESMFold 普遍 40-50, 用 35-55 区间映射
    plddt_score = max(0, min(1, (plddt - 35) / 20))

    # 2. chromophore 区域 pLDDT (更直接反映折叠正确性)
    #    cb 30-65 区间映射
    chromo_score = max(0, min(1, (chromo - 30) / 35))

    # 3. pTM 全局拓扑 (0.3-0.6 区间映射)
    ptm_score = max(0, min(1, (ptm - 0.3) / 0.3))

    # 4. 预期 Tm (热稳奖关键, Best Top-1 因子)
    #    Tm 70-95°C 映射, 72°C 比赛热处理温度
    tm_score = max(0, min(1, (tm - 70) / 25))

    # 5. 亮度风险 (突变数越多, Finit < 0.3WT 风险越高)
    #    n_mut 0-10 反向映射, mut>10 急剧惩罚
    if n_mut <= 5:
        brightness_score = 1.0
    elif n_mut <= 10:
        brightness_score = 1.0 - (n_mut - 5) * 0.05
    else:
        brightness_score = max(0, 0.75 - (n_mut - 10) * 0.15)

    # 综合权重 (按 Best Top-1 视角)
    weights = {
        "plddt": 1.5,        # 折叠是基础
        "chromo": 2.0,       # chromophore 区域更关键
        "ptm": 1.0,
        "tm": 2.5,           # 热稳奖直接乘法因子, 比赛核心
        "brightness": 2.0,   # 亮度阈值
    }

    total = (
        weights["plddt"] * plddt_score +
        weights["chromo"] * chromo_score +
        weights["ptm"] * ptm_score +
        weights["tm"] * tm_score +
        weights["brightness"] * brightness_score
    )
    max_total = sum(weights.values())  # 9.0
    score_normalized = total / max_total * 10  # 0-10

    return {
        "score": round(score_normalized, 2),
        "components": {
            "plddt_score": round(plddt_score, 3),
            "chromo_score": round(chromo_score, 3),
            "ptm_score": round(ptm_score, 3),
            "tm_score": round(tm_score, 3),
            "brightness_score": round(brightness_score, 3),
        }
    }


# ============================================================
# 计算所有得分
# ============================================================
print("=" * 105)
print(f"{'name':<32} {'role':<22} {'pLDDT':>5} {'cb':>5} {'pTM':>5} {'Tm':>3} {'mut':>3} {'score':>5}")
print("=" * 105)

for r in results:
    s = score_candidate(r)
    r["round4_score"] = s["score"]
    r["score_components"] = s["components"]

# 排序
sorted_results = sorted(results, key=lambda x: -x["round4_score"])
for r in sorted_results:
    print(f"{r['name']:<32} {r['role'][:22]:<22} "
          f"{r['plddt_mean']:>5.1f} {r['plddt_chromo_region']:>5.1f} "
          f"{r['ptm'] or 0:>5.3f} {r['expected_tm']:>3} {r['n_muts']:>3} {r['round4_score']:>5.2f}")

# ============================================================
# Top-6 选择策略 (角色多样化 + Best Top-1)
# ============================================================
print("\n" + "=" * 70)
print("Top-6 选择策略 (角色多样化 + Best Top-1):")
print("=" * 70)

# 按角色分组取最佳
by_role = defaultdict(list)
for r in sorted_results:
    by_role[r["role"]].append(r)

# 角色配额 (Best Top-1 视角):
quota = {
    "thermostable_hero": 1,    # mBaoJin 1 条 (上限高, 风险大)
    "combined_balanced": 2,    # htFuncLib 2 条 (主力)
    "safety_baseline": 1,      # sfGFP+I152S 1 条 (保险)
    "exploration": 2,          # 2 条探索
}

selected = []
for role, n in quota.items():
    role_picks = by_role[role][:n]
    selected.extend(role_picks)
    print(f"\n[{role}] 配额 {n}:")
    for r in role_picks:
        print(f"  ✓ {r['name']:<32} score={r['round4_score']:.2f}  Tm={r['expected_tm']}  mut={r['n_muts']}")

print(f"\n总选: {len(selected)} 条")

# ============================================================
# 优化: 如果某些条相同/相似 (例如 B2/B3 序列相同), 重新选择
# ============================================================
# 检查序列唯一性
seen_seqs = set()
deduped = []
for r in selected:
    if r["seq"] not in seen_seqs:
        seen_seqs.add(r["seq"])
        deduped.append(r)
    else:
        print(f"⚠ 序列重复, 跳过: {r['name']}")

# 补足 6 条
if len(deduped) < 6:
    print(f"\n选中去重后 {len(deduped)} 条, 需补足到 6...")
    remaining_pool = [r for r in sorted_results if r not in deduped and r["seq"] not in seen_seqs]
    for r in remaining_pool:
        if len(deduped) >= 6:
            break
        deduped.append(r)
        seen_seqs.add(r["seq"])
        print(f"  ➕ 补 {r['name']:<32} score={r['round4_score']:.2f}")

# 截断到 6
final_6 = deduped[:6]

# ============================================================
# 最终 6 条按 score 重排 (Seq_ID = 1 应是最强的)
# ============================================================
final_6_sorted = sorted(final_6, key=lambda x: -x["round4_score"])

print("\n" + "=" * 90)
print("Round 4 最终 6 条 (按综合得分降序):")
print("=" * 90)
print(f"{'Seq_ID':<6} {'name':<32} {'role':<22} {'score':>5} {'Tm':>3} {'mut':>3}")
print("-" * 90)
for i, r in enumerate(final_6_sorted, 1):
    print(f"{i:<6} {r['name']:<32} {r['role'][:22]:<22} {r['round4_score']:>5.2f} {r['expected_tm']:>3} {r['n_muts']:>3}")

# ============================================================
# 保存最终选择
# ============================================================
output_data = []
for i, r in enumerate(final_6_sorted, 1):
    output_data.append({
        "Seq_ID": i,
        "name": r["name"],
        "role": r["role"],
        "scaffold": r["scaffold"],
        "n_muts": r["n_muts"],
        "expected_tm": r["expected_tm"],
        "plddt_mean": r["plddt_mean"],
        "plddt_chromo_region": r["plddt_chromo_region"],
        "ptm": r["ptm"],
        "round4_score": r["round4_score"],
        "score_components": r["score_components"],
        "notes": r["notes"],
        "seq": r["seq"],
    })

with open(OUT / "final_6_round4.json", "w", encoding="utf-8") as f:
    json.dump(output_data, f, indent=2, ensure_ascii=False)

# 生成提交 CSV
sub = pd.DataFrame({
    "Team_Name": ["YourTeamName"] * 6,
    "Seq_ID": list(range(1, 7)),
    "Sequence": [r["seq"] for r in final_6_sorted]
})
sub_path = OUT / "submission_round4.csv"
sub.to_csv(sub_path, index=False)

print(f"\n✓ 最终 6 条保存到 work/round4/final_6_round4.json")
print(f"✓ 提交 CSV 保存到 {sub_path}")

# ============================================================
# 提交合规性最终检查
# ============================================================
print("\n" + "=" * 70)
print("提交合规性最终检查:")
print("=" * 70)
excl_seqs = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())
print(f"  Exclusion_List 共 {len(excl_seqs)} 条")
all_ok = True
for _, row in sub.iterrows():
    s = row["Sequence"]
    issues = []
    if not s.startswith("M"): issues.append("not-M-start")
    if not (220 <= len(s) <= 250): issues.append(f"len={len(s)}")
    if set(s) - set("ACDEFGHIKLMNPQRSTVWY"): issues.append("bad-AA")
    if "*" in s: issues.append("stop-codon")
    if s in excl_seqs: issues.append("EXCLUDED!")
    status = "✓ OK" if not issues else "❌ " + ", ".join(issues)
    print(f"  Seq {row['Seq_ID']}: len={len(s)}  {status}")
    if issues: all_ok = False

if all_ok:
    print("\n🎉 全部 6 条通过最终验证, 可提交!")
else:
    print("\n⚠️ 有问题需要修复!")
