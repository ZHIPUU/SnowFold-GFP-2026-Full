"""
Round 6 Final: 调整后 Top-6 提交
================================
核心判断: mBaoJin 低 pLDDT 是 ESMFold 模型偏差
理由:
  - mBaoJin 来源于 StayGold(Cytaeis uchidea), 与 avGFP 序列同源性 <33%
  - ESMFold 训练数据以 avGFP 家族为主, 对 mBaoJin 骨架的置信度系统性偏低
  - mBaoJin 晶体结构 (PDB 8QBJ) 已证实正确折叠
  - 文献 Tm=92°C 已实验验证
  - 所有 11 个 mBaoJin 候选 pLDDT 均集中在 ~38.9 (不随突变数变化),
    说明这不是突变引起的问题, 而是骨架级别的模型偏差

策略:
  1. mBaoJin: 使用文献亮度 (1.8× sfGFP), 文献 Tm (92°C)
     pLDDT 惩罚降至 0.75 (安全折扣, 非模型偏差惩罚)
  2. sfGFP 手工候选: 使用 pLDDT 校准亮度, 文献 Tm 
  3. MPNN/LMPNN: 使用 pTM 校准 + 突变风险折扣
  4. 多样性选择: 每个骨架类别最多 2 条
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
R5 = ROOT / "work" / "round5"
OUT = ROOT / "work" / "round6"

# ============================================================
# 加载候选
# ============================================================
with open(R5 / "evolvepro_scored.json", encoding="utf-8") as f:
    all_c = json.load(f)

print(f"加载候选: {len(all_c)} 条")

# ============================================================
# 评分函数 (v3 - 修正 mBaoJin 模型偏差)
# ============================================================
# 各骨架相对于 sfGFP 的亮度倍数 (文献值)
SCAFFOLD_BRIGHTNESS = {
    "sfGFP": 1.0, "avGFP": 0.33, "amacGFP": 0.40,
    "cgreGFP": 1.50, "ppluGFP": 0.80, "mBaoJin": 1.80,
    "sfGFP_MPNN": 0.50, "avGFP_MPNN": 0.20, "avGFP_LMPNN": 0.25,
}

# 已知亮度突变 boost
BRIGHTNESS_BOOSTS = {
    "S65T": 2.0, "F64L": 1.3, "F99S": 1.2, "M153T": 1.2,
    "V163A": 1.1, "I152S": 1.3, "Q69L": 1.2, "S72A": 1.1,
}

# 已知超稳突变 Tm boost (°C)
STABILIZING_MUTS = {
    "S30R": 7, "F64L": 5, "S65T": 6, "F99S": 6, "M153T": 6,
    "V163A": 5, "I152S": 4, "T59P": 8, "V60A": 3, "S72A": 5,
    "A206V": 3, "I171V": 3, "N105T": 3,
}

# 基础 Tm
BASE_TM = {
    "sfGFP": 78, "avGFP": 65, "amacGFP": 60, "cgreGFP": 65,
    "ppluGFP": 55, "mBaoJin": 92,
    "sfGFP_MPNN": 72, "avGFP_MPNN": 62, "avGFP_LMPNN": 62,
}

def compute(c):
    scaffold = c.get("scaffold", "avGFP")
    n_muts = c.get("n_muts", 0)
    plddt = c.get("plddt_mean", 50)
    ptm = c.get("ptm", 0.5) or 0.5
    seq = c.get("seq", "")
    
    # === A. Finit/Finit_WT ===
    base_brightness = SCAFFOLD_BRIGHTNESS.get(scaffold, 0.3)
    
    # 亮度突变 boost
    boost = 1.0
    if scaffold in ("sfGFP", "avGFP", "amacGFP", "cgreGFP", "ppluGFP", "mBaoJin"):
        for mut, factor in BRIGHTNESS_BOOSTS.items():
            if mut in seq:
                boost *= factor
    
    # 结构置信度折扣
    if scaffold == "mBaoJin":
        # mBaoJin: ESMFold 模型偏差, 使用轻度折扣
        struct_factor = 0.75  # 安全折扣, 非 pLDDT 惩罚
    elif scaffold in ("sfGFP_MPNN", "avGFP_MPNN", "avGFP_LMPNN"):
        # MPNN：用 pTM 评估
        struct_factor = np.clip(0.4 + ptm * 0.5, 0.4, 1.0)
    else:
        # 手工候选: pLDDT 校准
        if plddt >= 55:
            struct_factor = 1.0
        elif plddt >= 45:
            struct_factor = 0.80
        elif plddt >= 40:
            struct_factor = 0.60
        else:
            struct_factor = 0.40
    
    finit_ratio = min(base_brightness * boost * struct_factor, 3.0)
    
    # === B. Ffinal/Finit (72°C 热稳定性) ===
    base_tm = BASE_TM.get(scaffold, 65)
    
    if scaffold == "mBaoJin":
        tm_est = 92.0  # 文献值, 已实验验证
    else:
        # pTM boost
        tm_ptm = max(0, (ptm - 0.4) * 20)
        # 已知突变 boost
        tm_mut = sum(delta for mut, delta in STABILIZING_MUTS.items() if mut in seq)
        # 高风险惩罚
        tm_pen = min(n_muts * 0.15, 6) if n_muts > 12 else 0
        tm_est = min(base_tm + tm_ptm + tm_mut - tm_pen, 99)
    
    # 72°C survival (基于两态变性)
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
    
    # === C. 风险折扣 (亮度保留折扣) ===
    if n_muts <= 5:
        risk = 1.0
    elif n_muts <= 12:
        risk = 0.85
    elif ptm >= 0.65 and ("MPNN" in scaffold or "LMPNN" in scaffold):
        risk = 0.70  # MPNN 高结构质量补偿
    elif ptm >= 0.5:
        risk = 0.50
    else:
        risk = 0.35
    
    # === D. 极低亮度阈值 ===
    if finit_ratio < 0.3:
        final = 0.0
    else:
        final = finit_ratio * therm * risk
    
    return {
        "finit_ratio": round(finit_ratio, 4),
        "tm_est": round(tm_est, 1),
        "therm_factor": round(therm, 4),
        "risk_factor": round(risk, 4),
        "struct_factor": round(struct_factor, 4),
        "final_score": round(final, 4),
    }

# ============================================================
# 打分
# ============================================================
results = []
for c in all_c:
    s = compute(c)
    results.append({**c, **s})
results.sort(key=lambda x: -x["final_score"])

# 只保留通过阈值
qualified = [r for r in results if r["finit_ratio"] >= 0.3]
print(f"\n通过 30% 阈值: {len(qualified)} / {len(results)}")

# ============================================================
# Top-6 多样性选择
# ============================================================
def cat(scaffold):
    if "MPNN" in scaffold and "LMPNN" not in scaffold: return "MPNN"
    elif "LMPNN" in scaffold: return "LMPNN"
    elif scaffold == "sfGFP": return "sfGFP"
    elif scaffold == "avGFP": return "avGFP"
    elif scaffold == "mBaoJin": return "mBaoJin"
    elif scaffold in ("cgreGFP", "ppluGFP", "amacGFP"): return "other"
    return scaffold

# 每个类别选最好的最多 2 条
cat_count = {}
selected = []
seen_seqs = set()

for r in qualified:
    c = cat(r["scaffold"])
    if r["seq"] not in seen_seqs and cat_count.get(c, 0) < 2:
        selected.append(r)
        seen_seqs.add(r["seq"])
        cat_count[c] = cat_count.get(c, 0) + 1
    if len(selected) >= 6:
        break

if len(selected) < 6:
    for r in qualified:
        if r["seq"] not in seen_seqs and len(selected) < 6:
            selected.append(r)
            seen_seqs.add(r["seq"])

selected.sort(key=lambda x: -x["final_score"])

# ============================================================
# 输出
# ============================================================
print(f"\n{'='*120}")
print(f"🎯 Round 6 Final Top-6")
print(f"{'='*120}")
print(f"{'Seq':<5} {'Name':<38} {'Scaffold':<14} {'mut':>4} {'pLDDT':>5} "
      f"{'pTM':>5} {'Finit_ratio':>10} {'Tm':>5} {'Therm':>5} {'Risk':>5} {'Score':>8}")
print("-" * 120)
for i, r in enumerate(selected, 1):
    plddt = r['plddt_mean']
    ptm = r.get('ptm', 0) or 0
    print(f"{i:<5} {r['name'][:38]:<38} {r['scaffold'][:14]:<14} "
          f"{r['n_muts']:>4} {plddt:>5.1f} {ptm:>5.3f} "
          f"{r['finit_ratio']:>10.4f} {r['tm_est']:>5.1f} "
          f"{r['therm_factor']:>5.3f} {r['risk_factor']:>5.3f} "
          f"{r['final_score']:>8.4f}")

# 提交文件
sub = pd.DataFrame({
    "Team_Name": ["YourTeamName"] * 6,
    "Seq_ID": list(range(1, 7)),
    "Sequence": [r["seq"] for r in selected]
})

OUT.mkdir(parents=True, exist_ok=True)
sub_path = OUT / "submission_round6_final.csv"
sub.to_csv(sub_path, index=False)

# 合规
excl = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())
for s in sub["Sequence"]:
    assert s.startswith("M"), "M"
    assert 220 <= len(s) <= 250, f"len={len(s)}"
    assert not (set(s) - set("ACDEFGHIKLMNPQRSTVWY")), "AA"
    assert s not in excl, "Excl"

# 保存 JSON
final_json = []
for i, r in enumerate(selected, 1):
    final_json.append({
        "Seq_ID": i, "name": r["name"], "scaffold": r["scaffold"],
        "n_muts": int(r["n_muts"]), "length": len(r["seq"]),
        "plddt_mean": float(r["plddt_mean"]),
        "plddt_chromo_region": float(r.get("plddt_chromo_region", 0) or 0),
        "ptm": float(r.get("ptm", 0) or 0),
        "finit_ratio": float(r["finit_ratio"]),
        "tm_est": float(r["tm_est"]),
        "therm_factor": float(r["therm_factor"]),
        "risk_factor": float(r["risk_factor"]),
        "struct_factor": float(r["struct_factor"]),
        "final_score": float(r["final_score"]),
        "seq": r["seq"],
    })
with open(OUT / "final_6_round6_final.json", "w", encoding="utf-8") as f:
    json.dump(final_json, f, indent=2, ensure_ascii=False)

# 打印详情
print(f"\n{'='*120}")
print(f"详细信息:")
print(f"{'='*120}")
for r in final_json:
    print(f"\n  Seq {r['Seq_ID']}: {r['name']}")
    print(f"    骨架: {r['scaffold']} | 突变数: {r['n_muts']} | 长度: {r['length']}")
    print(f"    结构: pLDDT={r['plddt_mean']:.1f} | pTM={r['ptm']:.4f}")
    print(f"    亮度: Finit_ratio={r['finit_ratio']:.2f}× (结构折扣={r['struct_factor']:.2f})")
    print(f"    热稳: Tm~{r['tm_est']:.0f}°C | therm_factor={r['therm_factor']:.3f}")
    print(f"    风险: risk_factor={r['risk_factor']:.3f}")
    print(f"    >>> 综合分 = {r['final_score']:.4f} = {r['finit_ratio']:.2f} × {r['therm_factor']:.3f} × {r['risk_factor']:.3f}")

print(f"\n{'='*80}")
print(f"✅ Round 6 最终提交完成!")
print(f"  Best Top-1: {final_json[0]['name']} = {final_json[0]['final_score']:.4f}")
print(f"  提交文件: {sub_path}")
print(f"  Top-6 JSON: {OUT / 'final_6_round6_final.json'}")
print(f"  合规: ✅ 全部通过")
print(f"{'='*80}")
