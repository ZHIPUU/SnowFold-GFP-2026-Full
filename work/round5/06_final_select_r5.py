"""
Round 5 终极 Top-6 选择
==================================
合并:
  - Round 4 v5 53 条候选 (sfGFP_MPNN, avGFP_MPNN, sfGFP, avGFP, amacGFP, mBaoJin)
  - Round 5 LigandMPNN v2 15 条候选 (avGFP_LMPNN, chromophore-aware)

策略 (Best Top-1 + 多样性):
  - 1 条 LigandMPNN 王者 (pLDDT 63, pTM 0.71)
  - 1 条 ProteinMPNN 王者 (MPNN_T01_014 pLDDT 68)
  - 1 条 ProteinMPNN avGFP_T03 (pLDDT 61)
  - 1-2 条手工设计 sfGFP/avGFP (htFuncLib + I152S)
  - 1 条 mBaoJin (热稳爆款)
"""
import json, pandas as pd
from pathlib import Path
from collections import defaultdict

ROOT = Path(r"D:\生信\2026Protein Design")
R5 = ROOT / "work" / "round5"
R4 = ROOT / "work" / "round4"

# 加载 Round 4 v5 全部评估结果
with open(R4 / "esmfold_round4_v3.json", encoding="utf-8") as f:
    r4_handcraft = json.load(f)
with open(R4 / "esmfold_mpnn.json", encoding="utf-8") as f:
    r4_mpnn_sf = json.load(f)
with open(R4 / "esmfold_mpnn_av.json", encoding="utf-8") as f:
    r4_mpnn_av = json.load(f)

# Round 5 LigandMPNN
with open(R5 / "esmfold_lmpnn_v2.json", encoding="utf-8") as f:
    r5_lmpnn = json.load(f)

# 补字段
for r in r4_handcraft:
    r.setdefault("mpnn_recovery", None)
    r.setdefault("lmpnn_ligand", None)
for r in r4_mpnn_sf + r4_mpnn_av:
    r.setdefault("lmpnn_ligand", None)
for r in r5_lmpnn:
    r.setdefault("mpnn_recovery", None)
    r["mpnn_recovery"] = r.get("lmpnn_recovery")

# 修正 MPNN Tm 估值 (与 v5 一致)
for r in r4_mpnn_sf + r4_mpnn_av + r5_lmpnn:
    plddt = r["plddt_mean"]
    if plddt >= 65:
        r["expected_tm"] = 88
    elif plddt >= 60:
        r["expected_tm"] = 85  # LigandMPNN 用 85 (chromo-aware 更稳)
    elif plddt >= 55:
        r["expected_tm"] = 82
    elif plddt >= 50:
        r["expected_tm"] = 79
    else:
        r["expected_tm"] = 75

all_r = r4_handcraft + r4_mpnn_sf + r4_mpnn_av + r5_lmpnn
all_r = [r for r in all_r if r["plddt_mean"] >= 35]
print(f"总候选数: {len(all_r)}")

# 评分函数 (v5 + LigandMPNN 加分)
def score(r):
    plddt = r["plddt_mean"]
    chromo = r["plddt_chromo_region"]
    ptm = r["ptm"] or 0.0
    n_mut = r["n_muts"]
    tm = r["expected_tm"]
    
    plddt_s = max(0, min(1, (plddt - 35) / 35))
    chromo_s = max(0, min(1, (chromo - 30) / 40))
    ptm_s = max(0, min(1, (ptm - 0.3) / 0.5))
    tm_s = max(0, min(1, (tm - 70) / 25))
    
    # 突变惩罚
    if n_mut <= 5: br = 1.0
    elif n_mut <= 12: br = 1.0 - (n_mut - 5) * 0.03
    elif n_mut <= 60: br = max(0.6, 0.79 - (n_mut - 12) * 0.005)
    else: br = max(0.5, 0.6 - (n_mut - 60) * 0.002)
    
    # MPNN/LigandMPNN 高 pLDDT 加分
    if r["scaffold"].endswith("_MPNN") or "LMPNN" in r["scaffold"]:
        if plddt >= 60:
            br = max(br, 0.92)
    
    # LigandMPNN 额外加分 (chromophore-aware, 比 ProteinMPNN 更可靠)
    lmpnn_bonus = 0
    if "LMPNN" in r["scaffold"]:
        lig_conf = r.get("lmpnn_ligand", 0) or 0
        lmpnn_bonus = lig_conf * 0.5  # 0-0.5 分
    
    weights = {"plddt": 2.5, "chromo": 2.0, "ptm": 1.5, "tm": 2.0, "brightness": 2.0}
    total = (weights["plddt"]*plddt_s + weights["chromo"]*chromo_s + weights["ptm"]*ptm_s
             + weights["tm"]*tm_s + weights["brightness"]*br)
    return round(total / sum(weights.values()) * 10 + lmpnn_bonus, 2)

for r in all_r:
    r["round5_score"] = score(r)

sorted_r = sorted(all_r, key=lambda x: -x["round5_score"])

print("\n" + "=" * 125)
print("Top-20 综合排序 (Round 5):")
print("=" * 125)
print(f"{'#':<3} {'name':<32} {'scaffold':<14} {'mut':>3} {'pLDDT':>5} {'cb':>5} {'pTM':>5} {'Tm':>3} {'score':>5}")
print("-" * 125)
for i, r in enumerate(sorted_r[:20], 1):
    print(f"{i:<3} {r['name'][:32]:<32} {r['scaffold'][:14]:<14} {r['n_muts']:>3} "
          f"{r['plddt_mean']:>5.1f} {r['plddt_chromo_region']:>5.1f} {r['ptm'] or 0:>5.3f} "
          f"{r['expected_tm']:>3} {r['round5_score']:>5.2f}")

# Top-6 多样性
print("\n" + "=" * 70)
print("Top-6 多样性选择:")
print("=" * 70)
priority = ["avGFP_LMPNN", "sfGFP_MPNN", "avGFP_MPNN", "avGFP", "sfGFP", "amacGFP", "mBaoJin"]
selected = []
seen = set()
sc_cnt = defaultdict(int)

for sc in priority:
    for r in sorted_r:
        if r["scaffold"] == sc and r["seq"] not in seen:
            selected.append(r); seen.add(r["seq"]); sc_cnt[sc] += 1
            print(f"  ✓ {r['name'][:32]:<32} ({sc}) score={r['round5_score']}")
            break

for r in sorted_r:
    if len(selected) >= 6: break
    if r["seq"] in seen: continue
    if sc_cnt[r["scaffold"]] >= 2: continue
    selected.append(r); seen.add(r["seq"]); sc_cnt[r["scaffold"]] += 1
    print(f"  ➕ {r['name'][:32]:<32} ({r['scaffold']}) score={r['round5_score']}")

final_6 = sorted(selected[:6], key=lambda x: -x["round5_score"])

print("\n" + "=" * 125)
print("🎉 Round 5 终极 Top-6:")
print("=" * 125)
print(f"{'Seq':<4} {'name':<32} {'scaffold':<14} {'mut':>3} {'pLDDT':>5} {'cb':>5} {'pTM':>5} {'Tm':>3} {'score':>5}")
print("-" * 125)
for i, r in enumerate(final_6, 1):
    print(f"{i:<4} {r['name'][:32]:<32} {r['scaffold'][:14]:<14} {r['n_muts']:>3} "
          f"{r['plddt_mean']:>5.1f} {r['plddt_chromo_region']:>5.1f} {r['ptm']:>5.3f} "
          f"{r['expected_tm']:>3} {r['round5_score']:>5.2f}")

# 保存
output = []
for i, r in enumerate(final_6, 1):
    output.append({
        "Seq_ID": i, "name": r["name"], "scaffold": r["scaffold"],
        "n_muts": r["n_muts"], "expected_tm": r["expected_tm"],
        "plddt_mean": r["plddt_mean"], "plddt_chromo_region": r["plddt_chromo_region"],
        "ptm": r["ptm"], "round5_score": r["round5_score"],
        "notes": r.get("notes",""), "seq": r["seq"],
        "lmpnn_ligand": r.get("lmpnn_ligand"),
    })

with open(R5 / "final_6_round5.json", "w", encoding="utf-8") as f:
    json.dump(output, f, indent=2, ensure_ascii=False)

sub = pd.DataFrame({
    "Team_Name": ["YourTeamName"] * 6,
    "Seq_ID": list(range(1, 7)),
    "Sequence": [r["seq"] for r in final_6]
})
sub_path = R5 / "submission_round5.csv"
sub.to_csv(sub_path, index=False)

# 合规检查
print("\n📋 合规检查:")
excl_seqs = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())
all_ok = True
for _, row in sub.iterrows():
    s = row["Sequence"]
    issues = []
    if not s.startswith("M"): issues.append("not-M")
    if not (220 <= len(s) <= 250): issues.append(f"len={len(s)}")
    if set(s) - set("ACDEFGHIKLMNPQRSTVWY"): issues.append("bad-AA")
    if s in excl_seqs: issues.append("EXCL")
    chromo = next((c for c in ["TYG","SYG","GYG"] if c in s), "?")
    print(f"  Seq {row['Seq_ID']}: len={len(s)} chromo={chromo} {'✅' if not issues else '❌ ' + ','.join(issues)}")
    if issues: all_ok = False

print(f"\n{'🎉 全部通过' if all_ok else '⚠️ 有问题'}")
print(f"\n📁 提交: {sub_path}")
