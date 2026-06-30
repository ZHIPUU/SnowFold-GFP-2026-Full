"""
Round 4 Step 3B: 综合打分 + Top-6 最终选择 (扩展候选池)
==================================
策略 (Best Top-1):
  - 1 条 mBaoJin (高Tm爆款, 上限92°C)
  - 1 条 sfGFP+I152S 保险 (Round 3 验证)
  - 2-3 条 htFuncLib风格综合分
  - 1-2 条 探索/教科书
"""
import json, pandas as pd
from pathlib import Path
from collections import defaultdict

ROOT = Path(r"D:\生信\2026Protein Design")
OUT = ROOT / "work" / "round4"

with open(OUT / "esmfold_round4_extended.json", encoding="utf-8") as f:
    results = json.load(f)

def score_candidate(r):
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
    score_normalized = total / sum(weights.values()) * 10
    return round(score_normalized, 2), {
        "plddt": round(plddt_score, 3), "chromo": round(chromo_score, 3),
        "ptm": round(ptm_score, 3), "tm": round(tm_score, 3), "brightness": round(brightness_score, 3)
    }

# 打分
for r in results:
    score, comp = score_candidate(r)
    r["round4_score"] = score
    r["score_components"] = comp

sorted_results = sorted(results, key=lambda x: -x["round4_score"])

print("=" * 110)
print(f"{'#':<3} {'name':<32} {'role':<22} {'pLDDT':>5} {'cb':>5} {'pTM':>5} {'Tm':>3} {'mut':>3} {'score':>5}")
print("=" * 110)
for i, r in enumerate(sorted_results, 1):
    print(f"{i:<3} {r['name']:<32} {r['role'][:22]:<22} "
          f"{r['plddt_mean']:>5.1f} {r['plddt_chromo_region']:>5.1f} "
          f"{r['ptm'] or 0:>5.3f} {r['expected_tm']:>3} {r['n_muts']:>3} {r['round4_score']:>5.2f}")

# ============================================================
# Top-6 策略 (角色多样化 + Best Top-1 视角 + 去重)
# ============================================================
print("\n" + "=" * 70)
print("Top-6 选择 (Best Top-1 视角):")
print("=" * 70)

by_role = defaultdict(list)
for r in sorted_results:
    by_role[r["role"]].append(r)

# 优先级配额
quota = {
    "safety_baseline": 1,       # C1 sfGFP+I152S (Round 3最佳)
    "combined_balanced": 2,     # htFuncLib 综合分主力
    "exploration": 1,           # 教科书或新颖
    "thermostable_hero": 1,     # mBaoJin 高Tm爆款
    # 还差1条, 从未选中的最高分中补
}

selected = []
seen_seqs = set()

for role, n in quota.items():
    role_picks = []
    for r in by_role[role]:
        if r["seq"] not in seen_seqs and len(role_picks) < n:
            role_picks.append(r)
            seen_seqs.add(r["seq"])
    selected.extend(role_picks)
    print(f"\n[{role}] 配额 {n}:")
    for r in role_picks:
        print(f"  ✓ {r['name']:<32} score={r['round4_score']:.2f}  Tm={r['expected_tm']}  mut={r['n_muts']}")

# 补足到 6 (从剩下分数最高且不重复的选)
remaining = [r for r in sorted_results if r["seq"] not in seen_seqs]
print(f"\n剩余 {len(remaining)} 条未选, 补到 6:")
while len(selected) < 6 and remaining:
    r = remaining.pop(0)
    selected.append(r)
    seen_seqs.add(r["seq"])
    print(f"  ➕ {r['name']:<32} score={r['round4_score']:.2f}")

# ============================================================
# 重排 (Seq_ID 1 = 最高分)
# ============================================================
final_6 = sorted(selected[:6], key=lambda x: -x["round4_score"])

print("\n" + "=" * 110)
print("🎉 Round 4 最终 6 条 (按综合得分降序排列):")
print("=" * 110)
print(f"{'Seq_ID':<6} {'name':<32} {'role':<22} {'score':>5} {'pLDDT':>6} {'cb':>5} {'Tm':>3} {'mut':>3}")
print("-" * 110)
for i, r in enumerate(final_6, 1):
    print(f"{i:<6} {r['name']:<32} {r['role'][:22]:<22} {r['round4_score']:>5.2f} "
          f"{r['plddt_mean']:>6.1f} {r['plddt_chromo_region']:>5.1f} "
          f"{r['expected_tm']:>3} {r['n_muts']:>3}")

# ============================================================
# 保存
# ============================================================
output_data = []
for i, r in enumerate(final_6, 1):
    output_data.append({
        "Seq_ID": i, "name": r["name"], "role": r["role"],
        "scaffold": r["scaffold"], "n_muts": r["n_muts"],
        "expected_tm": r["expected_tm"],
        "plddt_mean": r["plddt_mean"], "plddt_chromo_region": r["plddt_chromo_region"],
        "ptm": r["ptm"], "round4_score": r["round4_score"],
        "score_components": r["score_components"],
        "notes": r["notes"], "seq": r["seq"],
    })

with open(OUT / "final_6_round4.json", "w", encoding="utf-8") as f:
    json.dump(output_data, f, indent=2, ensure_ascii=False)

sub = pd.DataFrame({
    "Team_Name": ["YourTeamName"] * 6,
    "Seq_ID": list(range(1, 7)),
    "Sequence": [r["seq"] for r in final_6]
})
sub_path = OUT / "submission_round4.csv"
sub.to_csv(sub_path, index=False)

print(f"\n✓ final_6_round4.json 保存")
print(f"✓ submission_round4.csv 保存到 {sub_path}")

# ============================================================
# 合规性检查
# ============================================================
print("\n" + "=" * 80)
print("📋 提交合规性最终检查 (官方规则):")
print("=" * 80)
excl_seqs = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())
all_ok = True
for _, row in sub.iterrows():
    s = row["Sequence"]
    issues = []
    if not s.startswith("M"): issues.append("not-M-start")
    if not (220 <= len(s) <= 250): issues.append(f"len={len(s)}_out_of_range")
    if set(s) - set("ACDEFGHIKLMNPQRSTVWY"): issues.append("bad-AA")
    if "*" in s: issues.append("stop-codon")
    if s in excl_seqs: issues.append("IN_EXCLUSION_LIST")
    status = "✅ OK" if not issues else "❌ " + ", ".join(issues)
    print(f"  Seq {row['Seq_ID']}: len={len(s)} starts={s[0]}  {status}")
    if issues: all_ok = False

if all_ok:
    print("\n🎉 全部 6 条满足官方提交要求！")
    print(f"\n📁 提交文件位置: {sub_path}")
    print("\n下一步: 写设计思路文档 + 建立开源仓库")
else:
    print("\n⚠️ 有合规性问题, 需要修复!")
