"""
Round 6 修正版: 结构完整性优先评分
====================================
核心原则 (来自用户):
  1. pTM < 0.5 = 结构崩塌 → 0 分 (mBaoJin 全系崩塌)
  2. pTM > 0.75 + 生色团 pLDDT > 80 = 金标结构 → 才可预测亮度
  3. MPNN 高 pTM 候选不应因突变数多而被惩罚
     → 它们就是被设计成这样的, 结构完整性才是价值所在
  4. 结构良好的 MPNN 至少能测到稳定基线荧光

当前现状:
  - 无候选通过金标 (pTM>0.75 + chromo>80)
  - 最接近: MPNN_T01_014 (pTM=0.765, chromo=67.6)
  - 这意味着还需更多轮设计来提升生色团质量

评分公式:
  综合分 = Finit/Finit_WT × Ffinal/Finit
  
  结构门控:
    pTM < 0.5 → 崩塌 → 0 分 (不进入评分)
    pTM >= 0.5 → 可能折叠 → 进入评分, 但用 pTM 校正亮度估计
  
  亮度估计:
    手工候选 (sfGFP/avGFP): 文献亮度 × pTM 校正
    MPNN/LMPNN: pLDDT 基准 × 生色团质量 × pTM 校正
    突变数不影响亮度估计 (MPNN 是整体设计的)
  
  热稳估计:
    pTM > 0.7 → 高折叠置信度 → 更高的 Tm 预期
    
提交策略:
  当前候选池中, 按结构质量排序选 Top-6
  最佳估计: 仍优于 Round 5 的手工打分
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
R5 = ROOT / "work" / "round5"
OUT = ROOT / "work" / "round6"

with open(R5 / "evolvepro_scored.json", encoding="utf-8") as f:
    all_c = json.load(f)

print(f"加载候选: {len(all_c)} 条")
print()

# ============================================================
# 正确的评分逻辑
# ============================================================

# 亮度基准 (文献)
KNOWN_BRIGHTNESS = {
    "sfGFP": 1.0, "avGFP": 0.33, "amacGFP": 0.40,
    "cgreGFP": 1.50, "ppluGFP": 0.80, "mBaoJin": 1.80,
    "sfGFP_MPNN": None, "avGFP_MPNN": None, "avGFP_LMPNN": None,  # 未知, 从结构推断
}

# 已知 Tm
KNOWN_TM = {
    "sfGFP": 78, "avGFP": 65, "amacGFP": 60, "cgreGFP": 65,
    "ppluGFP": 55, "mBaoJin": 92,
}

def score_candidate(c):
    scaffold = c.get("scaffold", "avGFP")
    n_muts = c.get("n_muts", 0)
    plddt = c.get("plddt_mean", 50) or 50
    chromo = c.get("plddt_chromo_region", plddt) or plddt
    ptm = c.get("ptm", 0.5) or 0.5
    seq = c.get("seq", "")
    
    # ============================================================
    # 第一关: 结构完整性门控
    # ============================================================
    if ptm < 0.5:
        return {
            "finit_ratio": 0.0, "tm_est": None,
            "therm_factor": 0.0,
            "plddt_grade": "COLLAPSED",
            "final_score": 0.0,
            "notes": f"pTM={ptm:.3f} < 0.5, 结构崩塌"
        }
    
    # 结构等级
    if ptm >= 0.75 and chromo >= 80:
        grade = "GOLD"
    elif ptm >= 0.7 and chromo >= 70:
        grade = "SILVER+"
    elif ptm >= 0.7:
        grade = "SILVER"
    elif ptm >= 0.65:
        grade = "BRONZE+"
    elif ptm >= 0.6:
        grade = "BRONZE"
    else:
        grade = "RISKY"
    
    # ============================================================
    # 第二关: Finit/Finit_WT (亮度估计)
    # ============================================================
    if scaffold in KNOWN_BRIGHTNESS and KNOWN_BRIGHTNESS[scaffold] is not None:
        # 已知骨架: 文献亮度 × 结构校正
        base = KNOWN_BRIGHTNESS[scaffold]
        # 结构校正: pTM 越高, 越接近文献亮度
        struct_confidence = 0.4 + 0.6 * np.clip((ptm - 0.5) / 0.3, 0, 1)
        # 生色团质量校正
        chromo_quality = np.clip(chromo / 80.0, 0.3, 1.0)
        finit = base * struct_confidence * chromo_quality
    else:
        # MPNN/LMPNN: 从结构推断
        # pLDDT 反映残基级置信度, 生色团 pLDDT 反映生色团质量
        if chromo >= 70:
            finit = 0.60  # 生色团区域好, 很可能有荧光
        elif chromo >= 60:
            finit = 0.40  # 生色团中等, 可能有弱荧光
        elif chromo >= 50:
            finit = 0.25  # 生色团一般
        else:
            finit = 0.15  # 生色团差
        
        # pTM 调整: pTM 高于 0.7 是加分
        ptm_bonus = 1.0 + max(0, (ptm - 0.7) * 1.5)
        finit *= ptm_bonus
    
    finit = np.clip(finit, 0.05, 3.0)
    
    # ============================================================
    # 30% 阈值检查
    # ============================================================
    if finit < 0.3:
        return {
            "finit_ratio": round(finit, 4), "tm_est": None,
            "therm_factor": 0.0,
            "plddt_grade": grade,
            "final_score": 0.0,
            "notes": f"Finit={finit:.2f}× < 30% WT → 淘汰"
        }
    
    # ============================================================
    # 第三关: Ffinal/Finit (72°C 热稳定性)
    # ============================================================
    if scaffold in KNOWN_TM:
        base_tm = KNOWN_TM[scaffold]
    else:
        base_tm = 70  # MPNN 默认
    
    # pTM 越高 → 结构越稳定 → Tm 越高
    tm_ptm = max(0, (ptm - 0.5) * 40)  # pTM=0.5→0, pTM=0.75→+10, pTM=0.8→+12
    tm_est = min(base_tm + tm_ptm, 99)
    
    # 72°C survival (两态模型简化)
    if tm_est >= 82:
        therm = 0.90 + 0.10 * (tm_est - 82) / 17
    elif tm_est >= 75:
        therm = 0.60 + 0.30 * (tm_est - 75) / 7
    elif tm_est >= 68:
        therm = 0.25 + 0.35 * (tm_est - 68) / 7
    elif tm_est >= 60:
        therm = 0.10 + 0.15 * (tm_est - 60) / 8
    else:
        therm = 0.03 + 0.07 * tm_est / 60
    therm = float(np.clip(therm, 0.03, 1.0))
    
    # ============================================================
    # 综合
    # ============================================================
    final_score = float(finit * therm)
    
    return {
        "finit_ratio": round(finit, 4),
        "tm_est": round(tm_est, 1),
        "therm_factor": round(therm, 4),
        "plddt_grade": grade,
        "final_score": round(final_score, 4),
        "notes": f"grade={grade}, ptm={ptm:.3f}, chromo_pLDDT={chromo:.1f}"
    }

# ============================================================
# 打分
# ============================================================
results = []
for c in all_c:
    s = score_candidate(c)
    results.append({**c, **s})

# 按综合分排序 (倒塌的直接 0 分在最后)
results.sort(key=lambda x: -x["final_score"])

print(f"{'='*130}")
print(f"修正后排名 (结构门控: pTM<0.5=0分)")
print(f"{'='*130}")
print(f"{'#':<4} {'Name':<35} {'Scaffold':<14} {'pLDDT':>5} {'chromo':>6} {'pTM':>6} "
      f"{'Grade':<10} {'Finit':>6} {'Tm':>5} {'Therm':>6} {'Score':>8}")
print("-" * 130)

gold_count = 0
for i, r in enumerate(results[:40], 1):
    ptm = r.get("ptm",0) or 0
    plddt = r.get("plddt_mean",0) or 0
    chromo = r.get("plddt_chromo_region", plddt) or plddt
    grade = r["plddt_grade"]
    if grade == "GOLD": gold_count += 1
    
    finit_display = r["finit_ratio"] if r["finit_ratio"] > 0 else 0
    score_display = r["final_score"] if r["final_score"] > 0 else 0
    therm_display = r["therm_factor"] if r["therm_factor"] > 0 else 0
    tm_display = r["tm_est"] if r["tm_est"] else 0
    
    print(f"{i:<4} {r['name'][:35]:<35} {r['scaffold'][:14]:<14} "
          f"{plddt:>5.1f} {chromo:>6.1f} {ptm:>6.4f} "
          f"{grade:<10} {finit_display:>6.3f} {tm_display:>5.1f} "
          f"{therm_display:>6.3f} {score_display:>8.4f}")

print(f"\n--- 金标 (GOLD) 候选: {gold_count} 条 ---")
print(f"注意: 当前无候选通过 pTM>0.75 + 生色团 pLDDT>80")
print(f"最接近: MPNN_T01_014 (pTM=0.765, chromo=67.6)")

# ============================================================
# Top-6: 按结构质量 + 亮度排, 多样性设计
# ============================================================
# 只保留结构合理的 (pTM >= 0.5 且 finit >= 0.3)
qualified = [r for r in results 
             if ((r.get("ptm",0) or 0) >= 0.5) and r["final_score"] > 0]
print(f"\n结构合理(pTM>=0.5)且通过阈值: {len(qualified)} / {len(results)}")

def cat(scaffold):
    if "MPNN" in scaffold and "LMPNN" not in scaffold: return "MPNN"
    elif "LMPNN" in scaffold: return "LMPNN"
    elif scaffold in ("sfGFP", "avGFP", "amacGFP", "cgreGFP", "ppluGFP"): return "manual"
    return scaffold

# 多样性选择: 覆盖 MPNN, LMPNN, manual 三种类别
selected, seen_seqs = [], set()
for cat_name in ["MPNN", "LMPNN", "manual"]:
    cat_count = 0
    for r in qualified:
        if cat(r["scaffold"]) == cat_name and r["seq"] not in seen_seqs and cat_count < 2:
            selected.append(r)
            seen_seqs.add(r["seq"])
            cat_count += 1
    if len(selected) >= 6:
        break

if len(selected) < 6:
    for r in qualified:
        if r["seq"] not in seen_seqs:
            selected.append(r)
            seen_seqs.add(r["seq"])
        if len(selected) >= 6:
            break

selected.sort(key=lambda x: -x["final_score"])

print(f"\n{'='*120}")
print(f"🎯 Round 6 Final Top-6 (修正版)")
print(f"{'='*120}")
print(f"{'Seq':<5} {'Name':<38} {'Scaffold':<14} {'pLDDT':>5} {'chromo':>6} "
      f"{'pTM':>6} {'Grade':<10} {'Finit':>6} {'Tm':>5} {'Score':>8}")
print("-" * 120)
for i, r in enumerate(selected, 1):
    ptm = r.get("ptm",0) or 0
    plddt = r.get("plddt_mean",0) or 0
    chromo = r.get("plddt_chromo_region", plddt) or plddt
    print(f"{i:<5} {r['name'][:38]:<38} {r['scaffold'][:14]:<14} "
          f"{plddt:>5.1f} {chromo:>6.1f} {ptm:>6.4f} "
          f"{r['plddt_grade']:<10} {r['finit_ratio']:>6.3f} "
          f"{r['tm_est'] or 0:>5.1f} {r['final_score']:>8.4f}")

# ============================================================
# 提交
# ============================================================
sub = pd.DataFrame({
    "Team_Name": ["YourTeamName"] * 6,
    "Seq_ID": list(range(1, 7)),
    "Sequence": [r["seq"] for r in selected]
})

OUT.mkdir(parents=True, exist_ok=True)
sub_path = OUT / "submission_round6_corrected.csv"
sub.to_csv(sub_path, index=False)

# 合规
excl = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())
for s in sub["Sequence"]:
    assert s.startswith("M"), "M"
    assert 220 <= len(s) <= 250, f"len={len(s)}"
    assert not (set(s) - set("ACDEFGHIKLMNPQRSTVWY")), "AA"
    assert s not in excl, "Excl"

print(f"\n{'='*80}")
print(f"✅ Round 6 修正版提交生成!")
print(f"  Best Top-1: {selected[0]['name']} = {selected[0]['final_score']:.4f}")
print(f"  提交文件: {sub_path}")
print(f"  合规: ✅ 全部通过")
print(f"{'='*80}")

# 与上一轮对比
print(f"\n{'='*80}")
print(f"对比 Round 5 v3 (原手工评分)")
print(f"{'='*80}")
with open(R5 / "final_6_round5_v3.json", encoding="utf-8") as f:
    r5_old = json.load(f)
print(f"  Round 5 v3 Top-1: {r5_old[0]['name']} pLDDT={r5_old[0]['plddt_mean']:.1f} pTM={r5_old[0]['ptm']:.4f}")
print(f"  Round 6 修正版 Top-1: {selected[0]['name']} pLDDT={selected[0].get('plddt_mean',0):.1f} pTM={selected[0].get('ptm',0):.4f}")
print(f"  核心改进: 结构门控淘汰了倒塌候选, MPNN 高 pTM 候选不被低估")
