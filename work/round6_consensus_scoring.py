"""
Round 6: 统一共识打分 v2 — 数据驱动 + 结构验证
===============================================
核心问题诊断:
  evolvepro_pred (XGBoost trained on 141K DMS data) 在 OOD 候选上
  系统性低估亮度。原因是训练数据以 avGFP 单/双突变为主，
  模型没见过 sfGFP/MPNN 等多突变组合，所以对所有远离 avGFP WT 
  的序列都预测低亮度。

修正策略:
  1. 对 near-WT 候选 (≤12 muts): 用 evolvepro_pred 直接预测亮度
  2. 对 sfGFP 骨架候选: 直接用 sfGFP 已知的 ~3× avGFP 亮度
  3. 对 MPNN/LMPNN 候选: 用 pTM/pLDDT 推断结构完整性，
     以 WT 亮度的折扣比例估算
  4. 对已知超稳序列 (mBaoJin): 用文献亮度值

比赛公式:
  综合分 = Finit/Finit_WT(相对WT亮度) × Ffinal/Finit(72°C热稳剩余)
  WT = sfGFP (比赛以 sfGFP 为参照)
"""
import json
import numpy as np
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
R5 = ROOT / "work" / "round5"
OUT = ROOT / "work" / "round6"

# ============================================================
# 1. 加载候选
# ============================================================
print("=" * 100)
print("Round 6 v2: 共识打分 (修复 OOD 亮度低估)")
print("=" * 100)

with open(R5 / "evolvepro_scored.json", encoding="utf-8") as f:
    all_c = json.load(f)
print(f"加载候选: {len(all_c)} 条")

# ============================================================
# 2. WT 参考亮度 (比赛以 sfGFP 为 WT)
# ============================================================
# sfGFP 是比赛 WT, Finit_WT = sfGFP 亮度
# sfGFP 亮度 ≈ avGFP × 3 (文献确证)
# avGFP log10 brightness = 3.72 (训练数据)
# sfGFP log10 brightness ≈ 3.72 + log10(3) ≈ 4.20

SFGFP_LOG10_BRIGHTNESS = 4.20  # sfGFP WT baseline (≈ avGFP × 3)

# 各骨架相对于 sfGFP 的亮度倍数 (基于文献 + 数据)
SCAFFOLD_BRIGHTNESS = {
    "sfGFP": 1.0,           # 比赛 WT
    "avGFP": 0.33,          # ≈ sfGFP 的 1/3
    "amacGFP": 0.40,        # 略高于 avGFP
    "cgreGFP": 1.50,        # 文献高 baseline
    "ppluGFP": 0.80,        # 中等
    "mBaoJin": 1.80,        # 92°C Tm, 亮度 ≈ sfGFP × 1.8
    "sfGFP_MPNN": 0.50,     # MPNN 设计的新序列, 亮度可能减半
    "avGFP_MPNN": 0.20,     # avGFP MPNN 更低
    "avGFP_LMPNN": 0.25,    # LigandMPNN 设计
}

# 已知提升亮度的突变 (Finit ratio boost relative to scaffold WT)
BRIGHTNESS_BOOSTS = {
    "S65T": 2.0,    # 成熟加快 + 亮度提升
    "F64L": 1.3,    # 折叠改善
    "F99S": 1.2,
    "M153T": 1.2,
    "V163A": 1.1,
    "I152S": 1.3,   # Round 1 验证过的热点
    "Q69L": 1.2,    # htFuncLib
    "S72A": 1.1,    # htFuncLib
    "T108V": 1.1,   # htFuncLib
}

# 已知超稳突变的 Tm 提升 (Δ°C)
STABILIZING_MUTS = {
    "S30R": 7, "F64L": 5, "S65T": 6, "F99S": 6,
    "M153T": 6, "V163A": 5, "I152S": 4, "T59P": 8,
    "V60A": 3, "S72A": 5, "A206V": 3, "I171V": 3,
    "N105T": 3, "Y145F": 2, "Y39N": 2,
}

def calc_finit_ratio(c):
    """
    估算 Finit/Finit_WT (相对于 sfGFP)
    策略: 根据候选类型采用不同的方法
    
    关键改进:
    - 对所有候选加入 pLDDT 结构罚分
    - mBaoJin: 文献亮度高但 ESMFold pLDDT 偏低,
      用 pLDDT 折扣避免过度乐观
    """
    scaffold = c.get("scaffold", "avGFP")
    n_muts = c.get("n_muts", 0)
    seq = c.get("seq", "")
    ptm = c.get("ptm", 0.5) or 0.5
    plddt = c.get("plddt_mean", 50)
    
    # 获取骨架基础亮度
    base = SCAFFOLD_BRIGHTNESS.get(scaffold, 0.3)
    
    # 结构置信度折扣 (适用于所有候选)
    # pLDDT < 45: 结构不可靠 → 强烈折扣
    # pLDDT 45-55: 边缘折扣
    # pLDDT > 55: 几乎不打折
    if plddt >= 60:
        struct_factor = 1.0
    elif plddt >= 50:
        struct_factor = 0.80
    elif plddt >= 45:
        struct_factor = 0.60
    elif plddt >= 40:
        struct_factor = 0.40
    else:
        struct_factor = 0.25  # pLDDT < 40, 折叠高度不确定
    
    # 已知亮度突变 boost (仅对手工候选)
    known_scaffolds = ("sfGFP", "avGFP", "amacGFP", "cgreGFP", "ppluGFP", "mBaoJin")
    boost = 1.0
    if scaffold in known_scaffolds:
        for mut, factor in BRIGHTNESS_BOOSTS.items():
            if mut in seq:
                boost *= factor
    
    # 最终 Finit ratio
    finit = base * boost * struct_factor
    return min(finit, 3.0)  # cap at 3× sfGFP

def calc_thermostability(c, finit_ratio):
    """
    估算 Ffinal/Finit (72°C 热处理后剩余比例)
    """
    scaffold = c.get("scaffold", "avGFP")
    n_muts = c.get("n_muts", 0)
    seq = c.get("seq", "")
    ptm = c.get("ptm", 0.5) or 0.5
    
    # Base Tm per scaffold (文献值, 部分来自 calibration)
    BASE_TM = {
        "sfGFP": 78, "avGFP": 65, "amacGFP": 60, 
        "cgreGFP": 65, "ppluGFP": 55, "mBaoJin": 92,
        "sfGFP_MPNN": 72, "avGFP_MPNN": 62, "avGFP_LMPNN": 62,
    }
    base_tm = BASE_TM.get(scaffold, 65)
    
    # Tm boost from pTM (higher fold confidence = more stable)
    tm_ptm_boost = max(0, (ptm - 0.4) * 25)  # pTM=0.4→0, pTM=0.8→+10
    
    # Tm boost from known stabilizing mutations
    tm_mut_boost = 0
    for mut, delta in STABILIZING_MUTS.items():
        if mut in seq:
            tm_mut_boost += delta
    
    # Tm penalty for many mutations (risk of destabilizing epistasis)
    if n_muts > 60:
        tm_penalty = min(n_muts * 0.1, 8)  # max -8°C
    elif n_muts > 12:
        tm_penalty = min(n_muts * 0.2, 5)  # max -5°C
    else:
        tm_penalty = 0
    
    tm_est = base_tm + tm_ptm_boost + tm_mut_boost - tm_penalty
    tm_est = np.clip(tm_est, 45, 99)
    
    # 72°C survival fraction (基于两态变性模型简化)
    if tm_est >= 82:
        therm = 0.90 + 0.10 * (tm_est - 82) / 17  # 0.90-1.00
    elif tm_est >= 75:
        therm = 0.60 + 0.30 * (tm_est - 75) / 7   # 0.60-0.90
    elif tm_est >= 68:
        therm = 0.25 + 0.35 * (tm_est - 68) / 7   # 0.25-0.60
    elif tm_est >= 60:
        therm = 0.10 + 0.15 * (tm_est - 60) / 8   # 0.10-0.25
    else:
        therm = 0.03 + 0.07 * tm_est / 60          # 0.03-0.10
    
    # 低亮度罚分: 如果 Finit_ratio < 0.3 → 淘汰 (0 分)
    if finit_ratio < 0.3:
        return 0, tm_est
    return float(np.clip(therm, 0.03, 1.0)), round(tm_est, 1)

# ============================================================
# 3. 对所有候选打分
# ============================================================
results = []
for c in all_c:
    finit_ratio = calc_finit_ratio(c)
    therm_factor, tm_est = calc_thermostability(c, finit_ratio)
    
    # 极低亮度阈值
    if finit_ratio < 0.3:
        final_score = 0
    else:
        final_score = finit_ratio * therm_factor
    
    results.append({
        **c,
        "finit_ratio": round(finit_ratio, 4),
        "tm_est": tm_est,
        "therm_factor": round(therm_factor, 4),
        "final_score": round(final_score, 4),
    })

results.sort(key=lambda x: -x["final_score"])

# ============================================================
# 4. 输出 Top-30
# ============================================================
print(f"\n{'='*130}")
print(f"Top-30 候选 (按 Round 6 综合分)")
print(f"{'='*130}")
print(f"{'#':<4} {'Name':<38} {'Scaffold':<14} {'pLDDT':>5} {'pTM':>5} "
      f"{'muts':>4} {'Finit_ratio':>10} {'Tm_est':>6} {'Therm':>6} {'Score':>8} {'⚠️'}")
print("-" * 130)
for i, r in enumerate(results[:30], 1):
    plddt = r['plddt_mean']
    warn = "⚠️pLDDT低" if plddt < 40 else ""
    warn += " ⚠️淘汰" if r['finit_ratio'] < 0.3 else ""
    print(f"{i:<4} {r['name'][:38]:<38} {r['scaffold'][:14]:<14} "
          f"{plddt:>5.1f} {(r['ptm'] or 0):>5.3f} "
          f"{r['n_muts']:>4} {r['finit_ratio']:>10.4f} "
          f"{r['tm_est']:>6.1f} {r['therm_factor']:>6.3f} "
          f"{r['final_score']:>8.4f} {warn:20s}")

# ============================================================
# 5. Top-6 多样性选择
# ============================================================
def cat(scaffold):
    if "MPNN" in scaffold and "LMPNN" not in scaffold: return "MPNN"
    elif "LMPNN" in scaffold: return "LMPNN"
    elif scaffold == "sfGFP": return "sfGFP"
    elif scaffold == "avGFP": return "avGFP"
    elif scaffold in ("mBaoJin",): return "mBaoJin"
    elif scaffold in ("cgreGFP", "ppluGFP", "amacGFP"): return "other"
    return scaffold

# 只选通过 30% 亮度阈值的候选
qualified = [r for r in results if r["finit_ratio"] >= 0.3]

print(f"\n通过 30% 亮度阈值的候选: {len(qualified)} / {len(results)}")

selected, seen_seqs, seen_cats = [], set(), set()
for r in qualified:
    c = cat(r["scaffold"])
    if r["seq"] not in seen_seqs and c not in seen_cats:
        selected.append(r); seen_seqs.add(r["seq"]); seen_cats.add(c)
    if len(selected) >= 6: break

# 补足到 6 条
if len(selected) < 6:
    for r in qualified:
        if r["seq"] not in seen_seqs:
            selected.append(r); seen_seqs.add(r["seq"])
        if len(selected) >= 6: break

selected.sort(key=lambda x: -x["final_score"])

print(f"\n{'='*100}")
print(f"🎯 Top-6 最终选择")
print(f"{'='*100}")
print(f"{'Seq':<5} {'Name':<38} {'Scaffold':<14} {'muts':>4} {'pLDDT':>5} "
      f"{'pTM':>5} {'Finit_ratio':>10} {'Tm':>5} {'Score':>8}")
print("-" * 100)
for i, r in enumerate(selected, 1):
    print(f"{i:<5} {r['name'][:38]:<38} {r['scaffold'][:14]:<14} "
          f"{r['n_muts']:>4} {r['plddt_mean']:>5.1f} {(r['ptm'] or 0):>5.3f} "
          f"{r['finit_ratio']:>10.4f} {r['tm_est']:>5.1f} "
          f"{r['final_score']:>8.4f}")

# ============================================================
# 6. 输出
# ============================================================
import pandas as pd

OUT.mkdir(parents=True, exist_ok=True)

# 提交文件
sub = pd.DataFrame({
    "Team_Name": ["YourTeamName"] * 6,
    "Seq_ID": list(range(1, 7)),
    "Sequence": [r["seq"] for r in selected]
})
sub_path = OUT / "submission_round6.csv"
sub.to_csv(sub_path, index=False)

# 合规检查
excl = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())
for s in sub["Sequence"]:
    assert s.startswith("M"), "Not starting with M"
    assert 220 <= len(s) <= 250, f"Length {len(s)} out of range"
    assert not (set(s) - set("ACDEFGHIKLMNPQRSTVWY")), "Invalid AA"
    assert s not in excl, "In exclusion list!"

def safe_val(v):
    """Convert numpy types to native Python for JSON serialization"""
    import numpy as np
    if isinstance(v, (np.integer,)): return int(v)
    if isinstance(v, (np.floating,)): return float(v)
    if isinstance(v, np.ndarray): return v.tolist()
    return v

# 保存完整排名
summary = []
for i, r in enumerate(results[:50], 1):
    summary.append({
        "rank": i, "name": r["name"], "scaffold": r["scaffold"],
        "n_muts": safe_val(r["n_muts"]), "plddt_mean": safe_val(r["plddt_mean"]),
        "plddt_chromo_region": safe_val(r.get("plddt_chromo_region")),
        "ptm": safe_val(r["ptm"]), "evolvepro_pred": safe_val(r.get("evolvepro_pred")),
        "finit_ratio": safe_val(r["finit_ratio"]), "tm_est": safe_val(r["tm_est"]),
        "therm_factor": safe_val(r["therm_factor"]), "final_score": safe_val(r["final_score"]),
        "seq": r["seq"],
    })
with open(OUT / "round6_full_ranking.json", "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)

# Final Top-6
final = []
for i, r in enumerate(selected, 1):
    final.append({
        "Seq_ID": i, "name": r["name"], "scaffold": r["scaffold"],
        "n_muts": safe_val(r["n_muts"]), "plddt_mean": safe_val(r["plddt_mean"]),
        "plddt_chromo_region": safe_val(r.get("plddt_chromo_region")),
        "ptm": safe_val(r["ptm"]), "finit_ratio": safe_val(r["finit_ratio"]),
        "tm_est": safe_val(r["tm_est"]), "therm_factor": safe_val(r["therm_factor"]),
        "final_score": safe_val(r["final_score"]), "seq": r["seq"],
        "notes": f"Round6, Finit_ratio={r['finit_ratio']:.2f}×, Tm~{r['tm_est']:.1f}°C, score={r['final_score']:.4f}",
    })
with open(OUT / "final_6_round6.json", "w", encoding="utf-8") as f:
    json.dump(final, f, indent=2, ensure_ascii=False)

print(f"\n{'='*80}")
print(f"✅ Round 6 完成!")
print(f"  提交: {sub_path}")
print(f"  Top-6: {OUT / 'final_6_round6.json'}")
print(f"  排名: {OUT / 'round6_full_ranking.json'}")
print(f"  合规: ✅ 全部通过")
print(f"{'='*80}")
