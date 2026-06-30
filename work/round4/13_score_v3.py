"""Round 4 v3 得分预测 (基于文献先验)"""
import json
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
OUT = ROOT / "work" / "round4"

with open(OUT / "final_6_round4_v3.json", encoding="utf-8") as f:
    final_6 = json.load(f)

def estimate_finit(c):
    sb = {"sfGFP": 1.00, "avGFP": 0.90, "amacGFP": 0.85, "cgreGFP": 1.10,
          "mBaoJin": 1.20, "sfGFP_MPNN": 1.00}
    base = sb.get(c["scaffold"], 1.0)
    
    mut_factor = 1.0
    name = c["name"]
    notes = c.get("notes", "")
    
    if c["scaffold"] == "avGFP" and ("sfGFP" in name or "sfGFP" in notes):
        mut_factor *= 1.5  # avGFP -> sfGFP 折叠改善大
    if "I152S" in name or "I152S" in notes:
        mut_factor *= 1.05
    if "Q69L" in name or "acid" in name.lower() or "Q69L" in notes:
        mut_factor *= 1.10
    if "S72A" in name:
        mut_factor *= 1.05
    if "S30R" in name:
        mut_factor *= 1.03
    if c["scaffold"] == "mBaoJin":
        mut_factor *= 0.98
    if c["scaffold"] == "sfGFP_MPNN":
        # MPNN de novo, 不确定性较大, 但 pLDDT 高代表折叠好
        mut_factor *= 1.0  # 中性
    
    plddt = c["plddt_mean"]
    if plddt >= 65: fold_prob = 0.98  # MPNN_T01_014 在此区间!
    elif plddt >= 55: fold_prob = 0.95
    elif plddt >= 45: fold_prob = 0.90
    elif plddt >= 40: fold_prob = 0.80
    elif plddt >= 35: fold_prob = 0.60
    else: fold_prob = 0.30
    
    chromo_p = c["plddt_chromo_region"]
    if chromo_p < 35: chromo_f = 0.70
    elif chromo_p < 45: chromo_f = 0.88
    elif chromo_p < 60: chromo_f = 0.96
    else: chromo_f = 1.00
    
    finit = base * mut_factor * fold_prob * chromo_f
    return finit, {"base": base, "mut_factor": mut_factor, "fold_prob": fold_prob, "chromo_factor": chromo_f}


def estimate_retention(c, test_T=72):
    tm = c["expected_tm"]
    delta = tm - test_T
    if delta >= 20: return 0.98
    if delta >= 14: return 0.92
    if delta >= 10: return 0.85
    if delta >= 6: return 0.70
    if delta >= 0: return 0.50
    return 0.20


print("=" * 105)
print(f"{'Seq':<4} {'name':<28} {'scaffold':<12} {'pLDDT':>5} {'cb':>5} {'Tm':>3} "
      f"{'Finit':>6} {'Therm':>6} {'综合':>5}")
print("=" * 105)

estimates = []
for c in final_6:
    finit, _ = estimate_finit(c)
    therm = estimate_retention(c)
    combined = finit * therm if finit >= 0.3 else 0
    icon = "🏆" if combined > 1.2 else ("⭐" if combined > 0.8 else "  ")
    print(f"{c['Seq_ID']:<4} {c['name'][:28]:<28} {c['scaffold'][:12]:<12} "
          f"{c['plddt_mean']:>5.1f} {c['plddt_chromo_region']:>5.1f} {c['expected_tm']:>3} "
          f"{finit:>6.2f} {therm:>6.2f} {combined:>5.2f} {icon}")
    estimates.append({"name": c["name"], "combined": combined, "finit": finit, "therm": therm})

top1 = max(estimates, key=lambda x: x["combined"])
print("\n" + "=" * 80)
print(f"🏆 Round 4 v3 预测 Best Top-1: {top1['name']}")
print(f"   综合分: {top1['combined']:.2f} (Finit {top1['finit']:.2f} × Therm {top1['therm']:.2f})")
print("=" * 80)

# 多场景
def pess(c):
    f, _ = estimate_finit(c); t = estimate_retention(c)
    return f * 0.65 * t * 0.85 if f >= 0.3 else 0
def opt(c):
    f, _ = estimate_finit(c); t = estimate_retention(c)
    return f * 1.4 * t * 1.10 if f >= 0.3 else 0

pess_s = [pess(c) for c in final_6]
neu_s = [e["combined"] for e in estimates]
opt_s = [opt(c) for c in final_6]

print(f"\n三场景 Best Top-1 预测:")
print(f"  🔴 悲观: {max(pess_s):.2f}")
print(f"  🟡 中性: {max(neu_s):.2f}")
print(f"  🟢 乐观: {max(opt_s):.2f}")

# 对比 v2 (前次预测)
print("\n" + "=" * 80)
print("📊 v2 vs v3 对比:")
print("=" * 80)
print(f"  v2 中性 Best Top-1: 1.36 (X4_avGFP_sfGFP_I152S_Q69L)")
print(f"  v3 中性 Best Top-1: {max(neu_s):.2f}")
delta = max(neu_s) - 1.36
print(f"  提升: {delta:+.2f} ({'+' if delta > 0 else ''}{delta/1.36*100:.0f}%)")
