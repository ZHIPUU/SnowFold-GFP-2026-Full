"""
Round 4 v5 - 修正 MPNN 的 Tm 估值并重新选择 Top-6
=====================================================
理论依据:
  - MPNN 设计旨在最大化结构-序列兼容性 (高 pLDDT)
  - 高 pLDDT 通常对应高 Tm (Sumida 2024 JACS, Dauparas 2022 Science)
  - sfGFP 同源蛋白 Tm 范围 ~78-96°C
  - MPNN 设计的 myoglobin 提升 Tm +27°C; TEV protease +20°C
  - 我们的 MPNN_T01_014 pLDDT 68.3 应对应较高 Tm

修正:
  - MPNN_T01_014 (sfGFP_MPNN, pLDDT 68): Tm 估 85°C (vs 之前 80)
  - MPNN_av_* (avGFP_MPNN, pLDDT 60+): Tm 估 82°C
  - 其他 MPNN < 50 pLDDT: 保持 78°C
"""
import json, pandas as pd
from pathlib import Path
from collections import defaultdict

ROOT = Path(r"D:\生信\2026Protein Design")
OUT = ROOT / "work" / "round4"

with open(OUT / "esmfold_round4_v3.json", encoding="utf-8") as f:
    handcraft = json.load(f)
with open(OUT / "esmfold_mpnn.json", encoding="utf-8") as f:
    mpnn_sf = json.load(f)
with open(OUT / "esmfold_mpnn_av.json", encoding="utf-8") as f:
    mpnn_av = json.load(f)

for r in handcraft:
    r.setdefault("mpnn_recovery", None); r.setdefault("mpnn_score", None)

# ============================================================
# 修正 MPNN Tm 估值 (基于 pLDDT)
# ============================================================
for r in mpnn_sf + mpnn_av:
    plddt = r["plddt_mean"]
    # MPNN 高 pLDDT -> 高 Tm
    if plddt >= 65:
        r["expected_tm"] = 88  # 类似 sfGFP+优化
    elif plddt >= 60:
        r["expected_tm"] = 84
    elif plddt >= 55:
        r["expected_tm"] = 81
    elif plddt >= 50:
        r["expected_tm"] = 78
    else:
        r["expected_tm"] = 72  # 风险大

all_r = [r for r in handcraft + mpnn_sf + mpnn_av if r["plddt_mean"] >= 35]
print(f"总候选: {len(all_r)}")

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
    else: br = 0.6
    if r["scaffold"].endswith("_MPNN") and plddt >= 60:
        br = max(br, 0.9)
    weights = {"plddt": 2.5, "chromo": 2.0, "ptm": 1.5, "tm": 2.0, "brightness": 2.0}
    total = (weights["plddt"]*plddt_s + weights["chromo"]*chromo_s + weights["ptm"]*ptm_s
             + weights["tm"]*tm_s + weights["brightness"]*br)
    return round(total / sum(weights.values()) * 10, 2)

for r in all_r:
    r["round4_score"] = score(r)

sorted_r = sorted(all_r, key=lambda x: -x["round4_score"])

print("\n" + "=" * 115)
print("Top-15 (修正 MPNN Tm 后):")
print("=" * 115)
print(f"{'#':<3} {'name':<28} {'scaffold':<12} {'mut':>3} {'pLDDT':>5} {'cb':>5} {'pTM':>5} {'Tm':>3} {'score':>5}")
print("-" * 115)
for i, r in enumerate(sorted_r[:15], 1):
    print(f"{i:<3} {r['name'][:28]:<28} {r['scaffold'][:12]:<12} {r['n_muts']:>3} "
          f"{r['plddt_mean']:>5.1f} {r['plddt_chromo_region']:>5.1f} {r['ptm'] or 0:>5.3f} "
          f"{r['expected_tm']:>3} {r['round4_score']:>5.2f}")

# Top-6 多样性
priority = ["sfGFP_MPNN", "avGFP_MPNN", "avGFP", "sfGFP", "amacGFP", "mBaoJin"]
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

final_6 = sorted(selected[:6], key=lambda x: -x["round4_score"])

print("\n" + "=" * 115)
print("🎉 Round 4 v5 终极版 Top-6 (修正Tm):")
print("=" * 115)
print(f"{'Seq':<4} {'name':<28} {'scaffold':<12} {'mut':>3} {'pLDDT':>5} {'cb':>5} {'pTM':>5} {'Tm':>3} {'score':>5}")
print("-" * 115)
for i, r in enumerate(final_6, 1):
    print(f"{i:<4} {r['name'][:28]:<28} {r['scaffold'][:12]:<12} {r['n_muts']:>3} "
          f"{r['plddt_mean']:>5.1f} {r['plddt_chromo_region']:>5.1f} {r['ptm']:>5.3f} "
          f"{r['expected_tm']:>3} {r['round4_score']:>5.2f}")

# ============================================================
# 比赛得分预测 (修正版)
# ============================================================
def estimate_finit(c):
    sb = {"sfGFP":1.00,"avGFP":0.90,"amacGFP":0.85,"cgreGFP":1.10,
          "mBaoJin":1.20,"sfGFP_MPNN":1.10,"avGFP_MPNN":1.00}
    base = sb.get(c["scaffold"], 1.0)
    mut_factor = 1.0
    name = c["name"]; notes = c.get("notes", "")
    if c["scaffold"] == "avGFP" and "sfGFP" in name: mut_factor *= 1.5
    if "I152S" in name: mut_factor *= 1.05
    if "Q69L" in name or "acid" in name.lower(): mut_factor *= 1.10
    if "S72A" in name: mut_factor *= 1.05
    if c["scaffold"] == "mBaoJin": mut_factor *= 0.98
    plddt = c["plddt_mean"]
    if plddt >= 65: fp = 0.98
    elif plddt >= 60: fp = 0.95
    elif plddt >= 55: fp = 0.92
    elif plddt >= 45: fp = 0.85
    elif plddt >= 40: fp = 0.75
    elif plddt >= 35: fp = 0.55
    else: fp = 0.25
    cb = c["plddt_chromo_region"]
    if cb < 35: cf = 0.65
    elif cb < 45: cf = 0.85
    elif cb < 55: cf = 0.95
    elif cb < 65: cf = 1.00
    else: cf = 1.05
    return base * mut_factor * fp * cf

def therm(c):
    d = c["expected_tm"] - 72
    if d >= 20: return 0.98
    if d >= 14: return 0.92
    if d >= 10: return 0.85
    if d >= 6: return 0.70
    if d >= 0: return 0.50
    return 0.20

print("\n" + "=" * 100)
print("💰 比赛得分预测 (基于修正 Tm):")
print("=" * 100)
print(f"{'Seq':<4} {'name':<28} {'Tm':>3} {'Finit':>6} {'Therm':>6} {'综合':>5}")
print("-" * 70)
estimates = []
for c in final_6:
    f = estimate_finit(c); t = therm(c)
    combined = f * t if f >= 0.3 else 0
    icon = "🏆" if combined > 1.2 else ("⭐" if combined > 0.9 else "  ")
    print(f"{c['Seq_ID'] if 'Seq_ID' in c else '-':<4} {c['name'][:28]:<28} "
          f"{c['expected_tm']:>3} {f:>6.2f} {t:>6.2f} {combined:>5.2f} {icon}")
    estimates.append(combined)

top1 = max(estimates)
print(f"\n🏆 v5 预测 Best Top-1: {top1:.2f}")

def pess(c):
    f = estimate_finit(c); t = therm(c)
    return f * 0.65 * t * 0.85 if f >= 0.3 else 0
def opt(c):
    f = estimate_finit(c); t = therm(c)
    return f * 1.4 * t * 1.10 if f >= 0.3 else 0

ps = [pess(c) for c in final_6]
os = [opt(c) for c in final_6]
print(f"\n三场景: 悲观 {max(ps):.2f} | 中性 {top1:.2f} | 乐观 {max(os):.2f}")

# 保存
output_data = []
for i, r in enumerate(final_6, 1):
    output_data.append({
        "Seq_ID": i, "name": r["name"], "scaffold": r["scaffold"],
        "n_muts": r["n_muts"], "expected_tm": r["expected_tm"],
        "plddt_mean": r["plddt_mean"], "plddt_chromo_region": r["plddt_chromo_region"],
        "ptm": r["ptm"], "round4_score": r["round4_score"],
        "notes": r.get("notes",""), "seq": r["seq"],
        "mpnn_recovery": r.get("mpnn_recovery"),
    })
with open(OUT / "final_6_round4_v5.json", "w", encoding="utf-8") as f:
    json.dump(output_data, f, indent=2, ensure_ascii=False)

sub = pd.DataFrame({
    "Team_Name": ["YourTeamName"] * 6,
    "Seq_ID": list(range(1, 7)),
    "Sequence": [r["seq"] for r in final_6]
})
sub_path = OUT / "submission_round4_v5.csv"
sub.to_csv(sub_path, index=False)
print(f"\n📁 提交: {sub_path}")

# 合规检查
excl_seqs = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())
print("\n📋 合规性:")
for _, row in sub.iterrows():
    s = row["Sequence"]
    ok = (s.startswith("M") and 220 <= len(s) <= 250 and
          not (set(s) - set("ACDEFGHIKLMNPQRSTVWY")) and s not in excl_seqs)
    chromo = next((c for c in ["TYG","SYG","GYG"] if c in s), "?")
    print(f"  Seq {row['Seq_ID']}: len={len(s)} chromo={chromo} {'✅' if ok else '❌'}")
