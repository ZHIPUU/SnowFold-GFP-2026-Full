"""
Round 5 紧急修正: 6 条低突变数安全候选 (全部 ≤13 mut)
==================================
策略转变:
  - 放弃 MPNN/LMPNN 高突变候选 (30% 阈值风险太大)
  - 用全部低突变手工设计 (突变 ≤13) 保证通过阈值
  - Best Top-1 从 sfGFP+htFuncLib 路线拿分
  
依据:
  - sfGFP WT Tm=86.1°C, 72°C 保留 ~85%
  - 突变 ≤5 的 sfGFP 变体亮度保留 >80% WT
  - htFuncLib sf:acid.3 论文验证 Tm 可达 96°C
"""
import json, pandas as pd
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
R5 = ROOT / "work" / "round5"

# 加载所有候选
with open(R5 / "evolvepro_scored.json", encoding="utf-8") as f:
    all_c = json.load(f)

# 只选突变数 ≤ 13 的 (安全通过 30% 阈值)
safe = [r for r in all_c if r["n_muts"] <= 13 and r["plddt_mean"] >= 35]
print(f"安全候选 (mut≤13, pLDDT≥35): {len(safe)}")

# 按预期综合分排序
# 综合 = 亮度比 × 热稳比
# 亮度比: mut≤5 → 0.8, mut 6-12 → 0.6, mut 13 → 0.5
# 热稳比: Tm 92 → 0.98, Tm 90 → 0.92, Tm 85 → 0.85, Tm 80 → 0.70
def estimate_combined(r):
    n_mut = r["n_muts"]
    tm = r.get("expected_tm", 80)
    
    # 修正 Tm (手工设计用文献值)
    plddt = r["plddt_mean"]
    if "MPNN" not in r["scaffold"] and "LMPNN" not in r["scaffold"]:
        # 手工设计用文献 Tm
        if r["scaffold"] == "sfGFP":
            if "acid" in r["name"] or "Q69L" in r["name"]:
                tm = 90  # htFuncLib 加成
            elif "S30R" in r["name"]:
                tm = 82
            else:
                tm = 80
        elif r["scaffold"] == "avGFP":
            if "sfGFP" in r["name"] and "I152S" in r["name"]:
                tm = 82  # Round 1 验证
            elif "sfGFP" in r["name"]:
                tm = 80
            else:
                tm = 75
        elif r["scaffold"] == "amacGFP":
            tm = 78
        elif r["scaffold"] == "mBaoJin":
            tm = 92
    
    # 亮度比
    if n_mut <= 5: finit = 0.85
    elif n_mut <= 10: finit = 0.65
    else: finit = 0.50
    
    # 热稳比
    delta = tm - 72
    if delta >= 20: therm = 0.98
    elif delta >= 14: therm = 0.92
    elif delta >= 10: therm = 0.85
    elif delta >= 6: therm = 0.70
    else: therm = 0.50
    
    return finit * therm, finit, therm, tm

for r in safe:
    combined, finit, therm, tm = estimate_combined(r)
    r["safe_combined"] = combined
    r["safe_finit"] = finit
    r["safe_therm"] = therm
    r["safe_tm"] = tm

safe.sort(key=lambda x: -x["safe_combined"])

print(f"\nTop-20 安全候选 (按预测综合分):")
print(f"{'#':<3} {'name':<35} {'scaffold':<14} {'mut':>3} {'pLDDT':>5} {'Tm':>3} {'Finit':>5} {'Therm':>5} {'综合':>5}")
print("-" * 110)
for i, r in enumerate(safe[:20], 1):
    print(f"{i:<3} {r['name'][:35]:<35} {r['scaffold'][:14]:<14} {r['n_muts']:>3} "
          f"{r['plddt_mean']:>5.1f} {r['safe_tm']:>3} {r['safe_finit']:>5.2f} "
          f"{r['safe_therm']:>5.2f} {r['safe_combined']:>5.2f}")

# Top-6 多样性选择
from collections import defaultdict
priority = ["sfGFP", "avGFP", "amacGFP", "mBaoJin", "cgreGFP"]
selected = []; seen = set(); sc_cnt = defaultdict(int)
for sc in priority:
    for r in safe:
        if r["scaffold"] == sc and r["seq"] not in seen:
            selected.append(r); seen.add(r["seq"]); sc_cnt[sc] += 1; break
for r in safe:
    if len(selected) >= 6: break
    if r["seq"] in seen: continue
    if sc_cnt[r["scaffold"]] >= 2: continue
    selected.append(r); seen.add(r["seq"]); sc_cnt[r["scaffold"]] += 1

final_6 = sorted(selected[:6], key=lambda x: -x["safe_combined"])

print(f"\n{'='*110}")
print(f"🛡️ Round 5 安全版 Top-6 (全部 mut≤13, 保证通过 30% 阈值):")
print(f"{'='*110}")
print(f"{'Seq':<4} {'name':<35} {'scaffold':<14} {'mut':>3} {'pLDDT':>5} {'Tm':>3} {'综合':>5}")
print("-" * 110)
for i, r in enumerate(final_6, 1):
    print(f"{i:<4} {r['name'][:35]:<35} {r['scaffold'][:14]:<14} {r['n_muts']:>3} "
          f"{r['plddt_mean']:>5.1f} {r['safe_tm']:>3} {r['safe_combined']:>5.2f}")

# 保存
output = []
for i, r in enumerate(final_6, 1):
    output.append({
        "Seq_ID": i, "name": r["name"], "scaffold": r["scaffold"],
        "n_muts": r["n_muts"], "expected_tm": r["safe_tm"],
        "plddt_mean": r["plddt_mean"], "plddt_chromo_region": r["plddt_chromo_region"],
        "ptm": r["ptm"], "safe_combined": r["safe_combined"],
        "safe_finit": r["safe_finit"], "safe_therm": r["safe_therm"],
        "evolvepro_pred": r.get("evolvepro_pred"),
        "notes": r.get("notes",""), "seq": r["seq"],
    })

with open(R5 / "final_6_round5_safe.json", "w", encoding="utf-8") as f:
    json.dump(output, f, indent=2, ensure_ascii=False)

sub = pd.DataFrame({
    "Team_Name": ["YourTeamName"] * 6,
    "Seq_ID": list(range(1, 7)),
    "Sequence": [r["seq"] for r in final_6]
})
sub_path = R5 / "submission_round5_safe.csv"
sub.to_csv(sub_path, index=False)

# 合规检查
excl = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())
all_ok = all(s.startswith("M") and 220<=len(s)<=250 and not (set(s)-set("ACDEFGHIKLMNPQRSTVWY")) and s not in excl for s in sub["Sequence"])
print(f"\n📋 合规: {'✅ 全部通过' if all_ok else '❌'}")
print(f"📁 提交: {sub_path}")

# Best Top-1
best = max(r["safe_combined"] for r in final_6)
print(f"\n🏆 Best Top-1 预测: {best:.2f}")
print(f"  (vs 当前 v3 的 0.34, 提升 {best - 0.34:+.2f})")
