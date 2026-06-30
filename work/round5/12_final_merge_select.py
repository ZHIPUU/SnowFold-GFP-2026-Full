"""Round 5 P0-3: 合并所有候选 + 最终 Top-6"""
import json, pandas as pd
from pathlib import Path
from collections import defaultdict

ROOT = Path(r"D:\生信\2026Protein Design")
R5 = ROOT / "work" / "round5"
R4 = ROOT / "work" / "round4"

# 加载所有
with open(R4 / "esmfold_round4_v3.json", encoding="utf-8") as f:
    r4 = json.load(f)
with open(R4 / "esmfold_mpnn.json", encoding="utf-8") as f:
    r4_sf = json.load(f)
with open(R4 / "esmfold_mpnn_av.json", encoding="utf-8") as f:
    r4_av = json.load(f)
with open(R5 / "esmfold_lmpnn_v2.json", encoding="utf-8") as f:
    r5_lm = json.load(f)
with open(R5 / "esmfold_lmpnn_expanded.json", encoding="utf-8") as f:
    r5_le = json.load(f)

# 补字段
for pool in [r4, r4_sf, r4_av, r5_lm, r5_le]:
    for r in pool:
        r.setdefault("mpnn_recovery", None)
        r.setdefault("lmpnn_ligand", None)
        r.setdefault("lmpnn_overall", None)

all_r = r4 + r4_sf + r4_av + r5_lm + r5_le
all_r = [r for r in all_r if r["plddt_mean"] >= 35]
print(f"总候选: {len(all_r)}")

# 修正 MPNN/LMPNN Tm
for r in all_r:
    plddt = r["plddt_mean"]
    if plddt >= 65: r["expected_tm"] = 88
    elif plddt >= 60: r["expected_tm"] = 85
    elif plddt >= 55: r["expected_tm"] = 82
    elif plddt >= 50: r["expected_tm"] = 79
    else: r["expected_tm"] = 75

# 评分
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
    lmpnn_bonus = 0
    if "LMPNN" in r["scaffold"]:
        lig = r.get("lmpnn_ligand", 0) or 0
        lmpnn_bonus = lig * 0.5
    w = {"plddt": 2.5, "chromo": 2.0, "ptm": 1.5, "tm": 2.0, "brightness": 2.0}
    total = (w["plddt"]*plddt_s + w["chromo"]*chromo_s + w["ptm"]*ptm_s + w["tm"]*tm_s + w["brightness"]*br)
    return round(total / sum(w.values()) * 10 + lmpnn_bonus, 2)

for r in all_r:
    r["round5_score"] = score(r)

sorted_r = sorted(all_r, key=lambda x: -x["round5_score"])

print("\nTop-15:")
print(f"{'#':<3} {'name':<35} {'scaffold':<14} {'mut':>3} {'pLDDT':>5} {'pTM':>5} {'Tm':>3} {'score':>5}")
print("-" * 100)
for i, r in enumerate(sorted_r[:15], 1):
    print(f"{i:<3} {r['name'][:35]:<35} {r['scaffold'][:14]:<14} {r['n_muts']:>3} "
          f"{r['plddt_mean']:>5.1f} {r['ptm'] or 0:>5.3f} {r['expected_tm']:>3} {r['round5_score']:>5.2f}")

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

final_6 = sorted(selected[:6], key=lambda x: -x["round5_score"])

print(f"\n🎉 Round 5 最终 Top-6:")
for i, r in enumerate(final_6, 1):
    print(f"  Seq {i}: {r['name'][:35]:<35} ({r['scaffold']}) pLDDT={r['plddt_mean']:.1f} Tm={r['expected_tm']} mut={r['n_muts']} score={r['round5_score']}")

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
with open(R5 / "final_6_round5_v2.json", "w", encoding="utf-8") as f:
    json.dump(output, f, indent=2, ensure_ascii=False)

sub = pd.DataFrame({
    "Team_Name": ["YourTeamName"] * 6,
    "Seq_ID": list(range(1, 7)),
    "Sequence": [r["seq"] for r in final_6]
})
sub_path = R5 / "submission_round5_v2.csv"
sub.to_csv(sub_path, index=False)

print(f"\n📁 提交: {sub_path}")
print("📋 合规: " + "✅ 全部通过" if all(220 <= len(r["seq"]) <= 250 and r["seq"].startswith("M") for r in final_6) else "❌")