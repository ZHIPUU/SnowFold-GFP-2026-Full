"""
Round 4 v4 终极版: 合并所有 65 条候选, 选 Top-6
==================================
候选池:
  - 41 手工设计 (sfGFP/avGFP/amacGFP/mBaoJin/cgreGFP)
  - 12 sfGFP_MPNN (de novo)
  - 8 avGFP_MPNN (de novo)
  
总计: 61 条 (剔除 cgreGFP 失败 + 4 重复)
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

# 补字段
for r in handcraft:
    r.setdefault("mpnn_recovery", None)
    r.setdefault("mpnn_score", None)

all_r = handcraft + mpnn_sf + mpnn_av
print(f"总评估: {len(all_r)}")

# 剔除明显失败
all_r = [r for r in all_r if r["plddt_mean"] >= 35]
print(f"过滤 pLDDT>=35: {len(all_r)}")

# 评分 (v4 调整: 拉高 pLDDT 权重)
def score(r):
    plddt = r["plddt_mean"]
    chromo = r["plddt_chromo_region"]
    ptm = r["ptm"] or 0.0
    n_mut = r["n_muts"]
    tm = r["expected_tm"]
    
    plddt_s = max(0, min(1, (plddt - 35) / 35))    # 35-70 区间, 高 pLDDT 拉满
    chromo_s = max(0, min(1, (chromo - 30) / 40))
    ptm_s = max(0, min(1, (ptm - 0.3) / 0.5))
    tm_s = max(0, min(1, (tm - 70) / 25))
    if n_mut <= 5: br = 1.0
    elif n_mut <= 12: br = 1.0 - (n_mut - 5) * 0.03
    elif n_mut <= 60: br = max(0.6, 0.79 - (n_mut - 12) * 0.005)  # MPNN 突变多但 pLDDT 弥补
    else: br = 0.6
    
    # MPNN 候选 pLDDT 高补偿
    if r["scaffold"].endswith("_MPNN"):
        if plddt >= 60:
            br = max(br, 0.9)  # pLDDT 充分证明折叠正确
    
    weights = {"plddt": 2.5, "chromo": 2.0, "ptm": 1.5, "tm": 2.0, "brightness": 2.0}
    total = (weights["plddt"]*plddt_s + weights["chromo"]*chromo_s + weights["ptm"]*ptm_s
             + weights["tm"]*tm_s + weights["brightness"]*br)
    return round(total / sum(weights.values()) * 10, 2)

for r in all_r:
    r["round4_score"] = score(r)

sorted_r = sorted(all_r, key=lambda x: -x["round4_score"])

print("\n" + "=" * 115)
print("Top-25 综合排序:")
print("=" * 115)
print(f"{'#':<3} {'name':<28} {'scaffold':<12} {'mut':>3} {'pLDDT':>5} {'cb':>5} {'pTM':>5} {'Tm':>3} {'score':>5}")
print("-" * 115)
for i, r in enumerate(sorted_r[:25], 1):
    print(f"{i:<3} {r['name'][:28]:<28} {r['scaffold'][:12]:<12} {r['n_muts']:>3} "
          f"{r['plddt_mean']:>5.1f} {r['plddt_chromo_region']:>5.1f} {r['ptm'] or 0:>5.3f} "
          f"{r['expected_tm']:>3} {r['round4_score']:>5.2f}")

# Top-6 多样性 (5 骨架优先, 每骨架最多 2)
print("\n" + "=" * 70)
print("Top-6 多样性选择:")
print("=" * 70)

priority = ["sfGFP_MPNN", "avGFP_MPNN", "avGFP", "sfGFP", "amacGFP", "mBaoJin"]
selected = []
seen = set()
sc_cnt = defaultdict(int)

# Phase 1: 每骨架最佳
for sc in priority:
    for r in sorted_r:
        if r["scaffold"] == sc and r["seq"] not in seen:
            selected.append(r); seen.add(r["seq"]); sc_cnt[sc] += 1
            print(f"  ✓ {r['name'][:28]:<28} ({sc}) score={r['round4_score']}")
            break

# Phase 2: 补到 6 (每骨架≤2)
for r in sorted_r:
    if len(selected) >= 6: break
    if r["seq"] in seen: continue
    if sc_cnt[r["scaffold"]] >= 2: continue
    selected.append(r); seen.add(r["seq"]); sc_cnt[r["scaffold"]] += 1
    print(f"  ➕ {r['name'][:28]:<28} ({r['scaffold']}) score={r['round4_score']}")

final_6 = sorted(selected[:6], key=lambda x: -x["round4_score"])

print("\n" + "=" * 115)
print("🎉 Round 4 v4 终极版 Top-6:")
print("=" * 115)
print(f"{'Seq':<4} {'name':<28} {'scaffold':<12} {'mut':>3} {'pLDDT':>5} {'cb':>5} {'pTM':>5} {'Tm':>3} {'score':>5}")
print("-" * 115)
for i, r in enumerate(final_6, 1):
    print(f"{i:<4} {r['name'][:28]:<28} {r['scaffold'][:12]:<12} {r['n_muts']:>3} "
          f"{r['plddt_mean']:>5.1f} {r['plddt_chromo_region']:>5.1f} {r['ptm']:>5.3f} "
          f"{r['expected_tm']:>3} {r['round4_score']:>5.2f}")

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
        "mpnn_score": r.get("mpnn_score"),
    })

with open(OUT / "final_6_round4_v4.json", "w", encoding="utf-8") as f:
    json.dump(output_data, f, indent=2, ensure_ascii=False)

sub = pd.DataFrame({
    "Team_Name": ["YourTeamName"] * 6,
    "Seq_ID": list(range(1, 7)),
    "Sequence": [r["seq"] for r in final_6]
})
sub_path = OUT / "submission_round4_v4.csv"
sub.to_csv(sub_path, index=False)

# 合规检查
print("\n📋 合规性检查:")
excl_seqs = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())
all_ok = True
for _, row in sub.iterrows():
    s = row["Sequence"]
    issues = []
    if not s.startswith("M"): issues.append("not-M")
    if not (220 <= len(s) <= 250): issues.append(f"len={len(s)}")
    if set(s) - set("ACDEFGHIKLMNPQRSTVWY"): issues.append("bad-AA")
    if s in excl_seqs: issues.append("EXCL!")
    chromo = next((c for c in ["TYG","SYG","GYG"] if c in s), "?")
    status = "✅" if not issues else "❌"
    print(f"  Seq {row['Seq_ID']}: len={len(s)} chromo={chromo} {status}" + (" " + ",".join(issues) if issues else ""))
    if issues: all_ok = False

print(f"\n{'🎉 全部通过' if all_ok else '⚠️ 有问题'}")
print(f"\n📁 提交: {sub_path}")
