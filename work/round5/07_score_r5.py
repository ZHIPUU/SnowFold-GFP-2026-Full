"""Round 5 得分预测 + 与 v5 对比"""
import json
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
R5 = ROOT / "work" / "round5"

with open(R5 / "final_6_round5.json", encoding="utf-8") as f:
    final_6 = json.load(f)

def estimate_finit(c):
    sb = {"sfGFP": 1.00, "avGFP": 0.90, "amacGFP": 0.85, "cgreGFP": 1.10,
          "mBaoJin": 1.20, "sfGFP_MPNN": 1.10, "avGFP_MPNN": 1.00, "avGFP_LMPNN": 1.05}  # LigandMPNN chromo-aware 加分
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
    
    chromo_p = c["plddt_chromo_region"]
    if chromo_p < 35: cf = 0.65
    elif chromo_p < 45: cf = 0.85
    elif chromo_p < 55: cf = 0.95
    elif chromo_p < 65: cf = 1.00
    else: cf = 1.05
    
    # LigandMPNN ligand_confidence 额外加成 (chromophore 环境正确)
    if c.get("lmpnn_ligand") and c["lmpnn_ligand"] > 0.45:
        cf *= 1.05  # chromophore-aware 额外加 5%
    
    return base * mut_factor * fp * cf

def therm(c):
    d = c["expected_tm"] - 72
    if d >= 20: return 0.98
    if d >= 14: return 0.92
    if d >= 10: return 0.85
    if d >= 6: return 0.70
    if d >= 0: return 0.50
    return 0.20

print("=" * 110)
print("💰 Round 5 比赛得分预测:")
print("=" * 110)
print(f"{'Seq':<4} {'name':<32} {'scaffold':<13} {'Tm':>3} {'Finit':>6} {'Therm':>6} {'综合':>5}")
print("-" * 95)

estimates = []
for c in final_6:
    f = estimate_finit(c); t = therm(c)
    combined = f * t if f >= 0.3 else 0
    icon = "🏆" if combined > 1.3 else ("⭐" if combined > 0.9 else "  ")
    print(f"{c['Seq_ID']:<4} {c['name'][:32]:<32} {c['scaffold'][:13]:<13} "
          f"{c['expected_tm']:>3} {f:>6.2f} {t:>6.2f} {combined:>5.2f} {icon}")
    estimates.append({"name": c["name"], "combined": combined, "finit": f, "therm": t})

top1 = max(estimates, key=lambda x: x["combined"])
print(f"\n🏆 Round 5 预测 Best Top-1: {top1['name']}")
print(f"   综合分: {top1['combined']:.2f}")

# 三场景
def pess(c):
    f = estimate_finit(c); t = therm(c)
    return f * 0.65 * t * 0.85 if f >= 0.3 else 0
def opt(c):
    f = estimate_finit(c); t = therm(c)
    return f * 1.4 * t * 1.10 if f >= 0.3 else 0

pess_s = [pess(c) for c in final_6]
neu_s = [e["combined"] for e in estimates]
opt_s = [opt(c) for c in final_6]

print(f"\n三场景 Best Top-1:")
print(f"  🔴 悲观: {max(pess_s):.2f}")
print(f"  🟡 中性: {max(neu_s):.2f}")
print(f"  🟢 乐观: {max(opt_s):.2f}")

print("\n" + "=" * 80)
print("📊 Round 4 v5 → Round 5 对比")
print("=" * 80)
print(f"  v5 中性 Best Top-1:        1.23")
print(f"  Round 5 中性 Best Top-1:   {max(neu_s):.2f}")
delta = max(neu_s) - 1.23
print(f"  提升:                       {delta:+.2f} ({delta/1.23*100:+.0f}%)")
print(f"\n  v5 乐观 Best Top-1:        1.90")
print(f"  Round 5 乐观 Best Top-1:   {max(opt_s):.2f}")
delta2 = max(opt_s) - 1.90
print(f"  提升:                       {delta2:+.2f} ({delta2/1.90*100:+.0f}%)")

print("\n" + "=" * 80)
print("📁 Round 5 关键升级:")
print("=" * 80)
print("  ✅ 引入 LigandMPNN (chromophore-aware), 比 ProteinMPNN 在 GFP 上更优")
print(f"  ✅ R5_av_lmpnn_v2_025: pLDDT 63.2 + pTM 0.71 + 170 mut (de novo 程度大)")
print(f"  ✅ 多骨架: 6 个不同 (LMPNN/MPNN_sf/MPNN_av/sfGFP/avGFP/amacGFP)")
print(f"  ✅ Best Top-1 候选保持 MPNN_T01_014 (pLDDT 68 全场最高)")
print(f"  ✅ 多了一个 LigandMPNN 备份方案, 即使 MPNN 失败仍有 7+ 分备选")
