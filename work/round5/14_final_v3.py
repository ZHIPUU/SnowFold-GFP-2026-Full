"""
Round 5 最终版: 综合 ESMFold pLDDT + EVOLVEpro + LigandMPNN 三维度评分
=================================================================
新信号 (Round 4 没有):
  1. EVOLVEpro: ESM-2 650M + XGBoost 顶层, val R²~0.5+ (与 Round 2 一致)
  2. LigandMPNN confidence: chromophore-aware 评分
  
关键调整:
  - EVOLVEpro 值低 (0.2-0.4) 但 *相对排序* 有意义
  - 候选间 evolvepro_pred 排名作为 brightness 排序参考
"""
import json, pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict

ROOT = Path(r"D:\生信\2026Protein Design")
R5 = ROOT / "work" / "round5"

# 加载完整候选 (含 evolvepro 评分)
with open(R5 / "evolvepro_scored.json", encoding="utf-8") as f:
    all_c = json.load(f)
print(f"候选 (含 EVOLVEpro 评分): {len(all_c)}")

# 过滤
all_c = [r for r in all_c if r["plddt_mean"] >= 35]
print(f"过滤 pLDDT>=35: {len(all_c)}")

# 修正 Tm
for r in all_c:
    plddt = r["plddt_mean"]
    if plddt >= 65: r["expected_tm"] = 88
    elif plddt >= 60: r["expected_tm"] = 85
    elif plddt >= 55: r["expected_tm"] = 82
    elif plddt >= 50: r["expected_tm"] = 79
    else: r["expected_tm"] = 75

# 计算 EVOLVEpro 百分位排名 (越高越好)
evo_scores = [r.get("evolvepro_pred", 0) for r in all_c]
evo_min = min(evo_scores); evo_max = max(evo_scores)
for r in all_c:
    r["evolvepro_normalized"] = (r.get("evolvepro_pred", 0) - evo_min) / (evo_max - evo_min + 1e-6)

# Round 5 综合评分
def score(r):
    plddt = r["plddt_mean"]; chromo = r["plddt_chromo_region"]
    ptm = r["ptm"] or 0.0; n_mut = r["n_muts"]; tm = r["expected_tm"]
    
    plddt_s = max(0, min(1, (plddt - 35) / 35))
    chromo_s = max(0, min(1, (chromo - 30) / 40))
    ptm_s = max(0, min(1, (ptm - 0.3) / 0.5))
    tm_s = max(0, min(1, (tm - 70) / 25))
    if n_mut <= 5: br = 1.0
    elif n_mut <= 12: br = 1.0 - (n_mut - 5) * 0.03
    elif n_mut <= 60: br = max(0.6, 0.79 - (n_mut - 12) * 0.005)
    else: br = max(0.5, 0.6 - (n_mut - 60) * 0.002)
    if r["scaffold"].endswith("_MPNN") or "LMPNN" in r["scaffold"]:
        if plddt >= 60: br = max(br, 0.92)
    
    # EVOLVEpro 加分 (新维度)
    evo_s = r["evolvepro_normalized"]
    
    # LigandMPNN ligand confidence (chromo-aware)
    lmpnn_bonus = 0
    if "LMPNN" in r["scaffold"]:
        lig = r.get("lmpnn_ligand", 0) or 0
        lmpnn_bonus = lig * 0.5
    
    w = {"plddt": 2.0, "chromo": 2.0, "ptm": 1.5, "tm": 2.0, "brightness": 2.0, "evo": 1.0}
    total = (w["plddt"]*plddt_s + w["chromo"]*chromo_s + w["ptm"]*ptm_s + w["tm"]*tm_s 
             + w["brightness"]*br + w["evo"]*evo_s)
    return round(total / sum(w.values()) * 10 + lmpnn_bonus, 2)

for r in all_c:
    r["round5_score_v3"] = score(r)

sorted_r = sorted(all_c, key=lambda x: -x["round5_score_v3"])

print("\nTop-20 (融合 ESMFold + EVOLVEpro + LigandMPNN):")
print(f"{'#':<3} {'name':<35} {'scaffold':<14} {'mut':>3} {'pLDDT':>5} {'pTM':>5} {'evo':>5} {'lig':>5} {'score':>5}")
print("-" * 110)
for i, r in enumerate(sorted_r[:20], 1):
    lig = r.get("lmpnn_ligand") or 0
    print(f"{i:<3} {r['name'][:35]:<35} {r['scaffold'][:14]:<14} {r['n_muts']:>3} "
          f"{r['plddt_mean']:>5.1f} {r['ptm'] or 0:>5.3f} {r['evolvepro_pred']:>5.2f} "
          f"{lig:>5.3f} {r['round5_score_v3']:>5.2f}")

# Top-6 多样性
priority = ["avGFP_LMPNN", "sfGFP_MPNN", "avGFP_MPNN", "avGFP", "sfGFP", "amacGFP", "mBaoJin"]
selected = []; seen = set(); sc_cnt = defaultdict(int)
for sc in priority:
    for r in sorted_r:
        if r["scaffold"] == sc and r["seq"] not in seen:
            selected.append(r); seen.add(r["seq"]); sc_cnt[sc] += 1; break
for r in sorted_r:
    if len(selected) >= 6: break
    if r["seq"] in seen: continue
    if sc_cnt[r["scaffold"]] >= 2: continue
    selected.append(r); seen.add(r["seq"]); sc_cnt[r["scaffold"]] += 1

final_6 = sorted(selected[:6], key=lambda x: -x["round5_score_v3"])

print(f"\n🎉 Round 5 v3 终极 Top-6:")
for i, r in enumerate(final_6, 1):
    lig = r.get("lmpnn_ligand") or 0
    print(f"  Seq {i}: {r['name'][:35]:<35} ({r['scaffold']:<12}) "
          f"pLDDT={r['plddt_mean']:.1f} pTM={r['ptm']:.3f} evo={r['evolvepro_pred']:.2f} "
          f"lig={lig:.3f} score={r['round5_score_v3']}")

# 保存
output = []
for i, r in enumerate(final_6, 1):
    output.append({
        "Seq_ID": i, "name": r["name"], "scaffold": r["scaffold"],
        "n_muts": r["n_muts"], "expected_tm": r["expected_tm"],
        "plddt_mean": r["plddt_mean"], "plddt_chromo_region": r["plddt_chromo_region"],
        "ptm": r["ptm"], "round5_score_v3": r["round5_score_v3"],
        "evolvepro_pred": r["evolvepro_pred"],
        "lmpnn_ligand": r.get("lmpnn_ligand"),
        "notes": r.get("notes",""), "seq": r["seq"],
    })
with open(R5 / "final_6_round5_v3.json", "w", encoding="utf-8") as f:
    json.dump(output, f, indent=2, ensure_ascii=False)

sub = pd.DataFrame({
    "Team_Name": ["YourTeamName"] * 6,
    "Seq_ID": list(range(1, 7)),
    "Sequence": [r["seq"] for r in final_6]
})
sub_path = R5 / "submission_round5_v3.csv"
sub.to_csv(sub_path, index=False)
print(f"\n📁 提交: {sub_path}")

# 合规检查
excl = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())
all_ok = all(s.startswith("M") and 220<=len(s)<=250 and not (set(s)-set("ACDEFGHIKLMNPQRSTVWY")) and s not in excl for s in sub["Sequence"])
print(f"📋 合规检查: {'✅ 全部通过' if all_ok else '❌ 失败'}")
