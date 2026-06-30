"""Round 4 v4 得分预测"""
import json
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
OUT = ROOT / "work" / "round4"

with open(OUT / "final_6_round4_v4.json", encoding="utf-8") as f:
    final_6 = json.load(f)

def estimate_finit(c):
    sb = {"sfGFP": 1.00, "avGFP": 0.90, "amacGFP": 0.85, "cgreGFP": 1.10,
          "mBaoJin": 1.20, "sfGFP_MPNN": 1.10, "avGFP_MPNN": 1.00}  # MPNN高pLDDT基础好
    base = sb.get(c["scaffold"], 1.0)
    
    mut_factor = 1.0
    name = c["name"]
    notes = c.get("notes", "")
    
    if c["scaffold"] == "avGFP" and ("sfGFP" in name or "sfGFP" in notes):
        mut_factor *= 1.5
    if "I152S" in name or "I152S" in notes:
        mut_factor *= 1.05
    if "Q69L" in name or "acid" in name.lower():
        mut_factor *= 1.10
    if "S72A" in name:
        mut_factor *= 1.05
    if c["scaffold"] == "mBaoJin":
        mut_factor *= 0.98
    
    plddt = c["plddt_mean"]
    if plddt >= 65: fold_prob = 0.98
    elif plddt >= 60: fold_prob = 0.95
    elif plddt >= 55: fold_prob = 0.92
    elif plddt >= 45: fold_prob = 0.85
    elif plddt >= 40: fold_prob = 0.75
    elif plddt >= 35: fold_prob = 0.55
    else: fold_prob = 0.25
    
    chromo_p = c["plddt_chromo_region"]
    if chromo_p < 35: chromo_f = 0.65
    elif chromo_p < 45: chromo_f = 0.85
    elif chromo_p < 55: chromo_f = 0.95
    elif chromo_p < 65: chromo_f = 1.00
    else: chromo_f = 1.05
    
    finit = base * mut_factor * fold_prob * chromo_f
    return finit

def therm(c):
    delta = c["expected_tm"] - 72
    if delta >= 20: return 0.98
    if delta >= 14: return 0.92
    if delta >= 10: return 0.85
    if delta >= 6: return 0.70
    if delta >= 0: return 0.50
    return 0.20

print("=" * 105)
print(f"{'Seq':<4} {'name':<28} {'scaffold':<12} {'pLDDT':>5} {'cb':>5} {'Tm':>3} {'Finit':>6} {'Therm':>6} {'综合':>5}")
print("=" * 105)

estimates = []
for c in final_6:
    f = estimate_finit(c); t = therm(c)
    combined = f * t if f >= 0.3 else 0
    icon = "🏆" if combined > 1.2 else ("⭐" if combined > 0.8 else "  ")
    print(f"{c['Seq_ID']:<4} {c['name'][:28]:<28} {c['scaffold'][:12]:<12} "
          f"{c['plddt_mean']:>5.1f} {c['plddt_chromo_region']:>5.1f} {c['expected_tm']:>3} "
          f"{f:>6.2f} {t:>6.2f} {combined:>5.2f} {icon}")
    estimates.append({"name": c["name"], "combined": combined, "finit": f, "therm": t})

top1 = max(estimates, key=lambda x: x["combined"])
print("\n" + "=" * 80)
print(f"🏆 Round 4 v4 预测 Best Top-1: {top1['name']}")
print(f"   综合分: {top1['combined']:.2f} (Finit {top1['finit']:.2f} × Therm {top1['therm']:.2f})")
print("=" * 80)

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
print("📊 v3 vs v4 对比:")
print("=" * 80)
print(f"  v3 中性 Best Top-1: 1.32 (H1_avGFP_sfGFP_acid3_I152S)")
print(f"  v4 中性 Best Top-1: {max(neu_s):.2f}")
delta = max(neu_s) - 1.32
pct = delta / 1.32 * 100
print(f"  提升: {delta:+.2f} ({'+' if delta > 0 else ''}{pct:.0f}%)")
print(f"\n  v3 乐观 Best Top-1: 2.03")
print(f"  v4 乐观 Best Top-1: {max(opt_s):.2f}")
delta2 = max(opt_s) - 2.03
pct2 = delta2 / 2.03 * 100
print(f"  提升: {delta2:+.2f} ({'+' if delta2 > 0 else ''}{pct2:.0f}%)")
