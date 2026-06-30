"""
Round 5 逐条候选实战风险评估
==================================
关键规则:
  - Finitial < 0.3 × Finitial_WT → 0 分 (淘汰)
  - 综合 = Finitial/Finitial_WT × Ffinal/Finitial
  - Best Top-1 排名

风险维度:
  1. 折叠概率 (pLDDT + 突变数)
  2. 亮度保持 (突变对 chromophore 影响)
  3. 热稳定性 (Tm vs 72°C)
  4. CFPS 表达兼容性
"""
import json
from pathlib import Path

with open(r"D:\生信\2026Protein Design\work\round5\final_6_round5_v3.json", encoding="utf-8") as f:
    final_6 = json.load(f)

print("=" * 120)
print("逐条候选实战风险评估")
print("=" * 120)

for c in final_6:
    print(f"\n{'='*100}")
    print(f"Seq {c['Seq_ID']}: {c['name']}")
    print(f"  骨架: {c['scaffold']}, 突变数: {c['n_muts']}, pLDDT: {c['plddt_mean']}, pTM: {c['ptm']}")
    print(f"  EVOLVEpro 预测亮度: {c['evolvepro_pred']:.3f} (log10)")
    print(f"{'='*100}")
    
    # 风险评估
    risks = []
    score_estimates = []
    
    # === 1. 折叠概率 ===
    plddt = c["plddt_mean"]
    n_mut = c["n_muts"]
    
    if plddt >= 65:
        fold_prob = 0.90
        risks.append("✅ pLDDT>65: 折叠概率高 (90%)")
    elif plddt >= 55:
        fold_prob = 0.75
        risks.append("🟡 pLDDT 55-65: 折叠概率中等 (75%)")
    elif plddt >= 45:
        fold_prob = 0.60
        risks.append("🟠 pLDDT 45-55: 折叠概率不确定 (60%)")
    else:
        fold_prob = 0.40
        risks.append("🔴 pLDDT<45: 折叠概率低 (40%)")
    
    # === 2. 突变数对亮度的影响 ===
    if n_mut <= 5:
        bright_risk = "极低"
        finit_ratio = 0.8  # 保留 80% 亮度
        risks.append(f"✅ 突变数 {n_mut}: 亮度保留极高 (>80% WT)")
    elif n_mut <= 12:
        bright_risk = "低"
        finit_ratio = 0.6
        risks.append(f"🟡 突变数 {n_mut}: 亮度保留中等 (~60% WT)")
    elif n_mut <= 60:
        bright_risk = "高"
        finit_ratio = 0.3
        risks.append(f"🟠 突变数 {n_mut}: 亮度可能大幅下降 (~30% WT)")
    else:
        bright_risk = "极高"
        finit_ratio = 0.15
        risks.append(f"🔴 突变数 {n_mut}: 亮度极可能 <30% WT (淘汰风险!)")
    
    # === 3. 30% 阈值风险 ===
    if finit_ratio < 0.3:
        threshold_risk = "极高"
        risks.append("🔴🔴 30% 阈值风险极高: 大概率被淘汰为 0 分")
    elif finit_ratio < 0.5:
        threshold_risk = "高"
        risks.append("🟠 30% 阈值风险高: 可能在边缘")
    elif finit_ratio < 0.8:
        threshold_risk = "中"
        risks.append("🟡 30% 阈值风险中等")
    else:
        threshold_risk = "低"
        risks.append("✅ 30% 阈值风险低")
    
    # === 4. 热稳定性 ===
    tm = c["expected_tm"]
    if tm >= 85:
        therm_retention = 0.92
        risks.append(f"✅ Tm~{tm}°C: 热稳定性保留 ~92%")
    elif tm >= 80:
        therm_retention = 0.85
        risks.append(f"🟡 Tm~{tm}°C: 热稳定性保留 ~85%")
    else:
        therm_retention = 0.70
        risks.append(f"🟠 Tm~{tm}°C: 热稳定性保留 ~70%")
    
    # === 综合预测 ===
    finit_rel = finit_ratio * fold_prob  # P(通过阈值) × 亮度比
    combined = finit_rel * therm_retention if finit_rel >= 0.3 else 0
    
    print(f"\n  风险评估:")
    for r in risks:
        print(f"    {r}")
    
    print(f"\n  📊 预测:")
    print(f"    折叠概率:     {fold_prob*100:.0f}%")
    print(f"    亮度保留比:   {finit_ratio*100:.0f}% WT")
    print(f"    通过30%阈值:  {'是' if finit_ratio >= 0.3 else '否 (0分!)'}")
    print(f"    热稳保留:     {therm_retention*100:.0f}%")
    print(f"    综合分预测:   {combined:.2f}" if combined > 0 else "    综合分预测:   0.00 (淘汰)")
    
    score_estimates.append(combined)

print(f"\n\n{'='*120}")
print("📊 6 条候选综合预测排名")
print("=" * 120)

# 重新计算
results = []
for c in final_6:
    plddt = c["plddt_mean"]
    n_mut = c["n_muts"]
    tm = c["expected_tm"]
    
    if plddt >= 65: fold_prob = 0.90
    elif plddt >= 55: fold_prob = 0.75
    elif plddt >= 45: fold_prob = 0.60
    else: fold_prob = 0.40
    
    if n_mut <= 5: finit_ratio = 0.8
    elif n_mut <= 12: finit_ratio = 0.6
    elif n_mut <= 60: finit_ratio = 0.3
    else: finit_ratio = 0.15
    
    if tm >= 85: therm = 0.92
    elif tm >= 80: therm = 0.85
    else: therm = 0.70
    
    finit_rel = finit_ratio * fold_prob
    combined = finit_rel * therm if finit_rel >= 0.3 else 0
    
    # 乐观场景 (MPNN 高 pLDDT 折叠更好)
    if "MPNN" in c["scaffold"] or "LMPNN" in c["scaffold"]:
        if plddt >= 60:
            fold_optimistic = 0.95
            finit_optimistic = max(finit_ratio, 0.4)  # 高 pLDDT 至少 40%
            combined_opt = finit_optimistic * fold_optimistic * therm
        else:
            combined_opt = combined * 1.3
    else:
        combined_opt = combined * 1.2
    
    results.append({
        "seq": c["Seq_ID"],
        "name": c["name"][:35],
        "scaffold": c["scaffold"],
        "n_mut": n_mut,
        "pLDDT": plddt,
        "neutral": combined,
        "optimistic": combined_opt,
        "pass_30pct": finit_ratio >= 0.3,
    })

results.sort(key=lambda x: -x["neutral"])

print(f"\n{'Seq':<4} {'name':<35} {'scaffold':<14} {'mut':>3} {'pLDDT':>5} {'中性':>6} {'乐观':>6} {'通过30%':>6}")
print("-" * 100)
for r in results:
    print(f"{r['seq']:<4} {r['name']:<35} {r['scaffold']:<14} {r['n_mut']:>3} "
          f"{r['pLDDT']:>5.1f} {r['neutral']:>6.2f} {r['optimistic']:>6.2f} "
          f"{'✅' if r['pass_30pct'] else '❌':>6}")

# Best Top-1
best_neutral = max(r["neutral"] for r in results)
best_optimistic = max(r["optimistic"] for r in results)
n_pass = sum(1 for r in results if r["pass_30pct"])

print(f"\n{'='*80}")
print(f"Best Top-1 预测:")
print(f"  中性场景: {best_neutral:.2f}")
print(f"  乐观场景: {best_optimistic:.2f}")
print(f"  通过 30% 阈值: {n_pass}/6 条")
print(f"{'='*80}")
