"""
严格对照新规则评估当前所有候选
====================================
规则硬门槛:
  1. 全局 pTM > 0.75 (竞争线 > 0.85)
  2. 全局 pLDDT > 80.0 (竞争线 > 85.0)
  3. 生色团核心(残基58~72及210~230) pLDDT > 85.0
  4. Initial/WT >= 0.30
  5. Final/Initial > 0.50
  
奖项线(综合分): 铜>0.30, 银>0.50, 金>0.80
"""
import json, sys
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
R5 = ROOT / "work" / "round5"

with open(R5 / "evolvepro_scored.json", encoding="utf-8") as f:
    all_c = json.load(f)

# sfGFP 的残基编号规则 (238 aa)
# 生色团核心区域: 58-72 (含 Y66-G67 生色团三联体) 和 210-230 (C端生色团结合袋)
CHROMO_REGIONS = [(58, 72), (210, 230)]

print("=" * 130)
print("严格评估所有候选 vs 规则硬门槛")
print("=" * 130)

# 逐个评估
grades = {"PASS_ALL": 0, "PASS_pTM": 0, "PASS_none": 0}
best = None
best_score = -1

for c in all_c:
    ptm = c.get("ptm", 0) or 0
    plddt = c.get("plddt_mean", 0) or 0
    chromo = c.get("plddt_chromo_region", plddt) or plddt
    name = c["name"]
    scaffold = c["scaffold"]
    
    # 规则检查
    pass_ptm = ptm > 0.75
    pass_plddt = plddt > 80.0
    pass_chromo = chromo > 85.0
    
    if pass_ptm: grades["PASS_pTM"] += 1
    
    # 模拟综合分 (用修正版方法估计)
    # 仅对通过结构门槛的候选估算
    if pass_ptm and pass_plddt and pass_chromo:
        grades["PASS_ALL"] += 1
    else:
        grades["PASS_none"] = grades.get("PASS_none", 0) + 1

    # 记录差距
    gap_ptm = max(0, 0.75 - ptm)
    gap_plddt = max(0, 80.0 - plddt)
    gap_chromo = max(0, 85.0 - chromo)
    total_gap = gap_ptm * 100 + gap_plddt + gap_chromo * 2
    
    # 排序用
    c["_total_gap"] = total_gap
    c["_gap_ptm"] = gap_ptm
    c["_gap_plddt"] = gap_plddt
    c["_gap_chromo"] = gap_chromo

# 按总差距排序 (最小差距优先)
all_c.sort(key=lambda x: x["_total_gap"])

print(f"\n{'Rank':<5} {'Name':<38} {'Scaffold':<14} {'pLDDT':>6} {'Chromo':>6} {'pTM':>6} "
      f"{'ΔpTM':>6} {'ΔpLDDT':>7} {'ΔChr':>7} {'Status':<20}")
print("-" * 130)

rank = 0
for c in all_c:
    ptm = c.get("ptm", 0) or 0
    plddt = c.get("plddt_mean", 0) or 0
    chromo = c.get("plddt_chromo_region", plddt) or plddt
    
    pass_ptm = ptm > 0.75
    pass_plddt = plddt > 80.0
    pass_chromo = chromo > 85.0
    
    if pass_ptm and pass_plddt and pass_chromo:
        status = "✅ 全部通过"
    elif pass_ptm:
        status = "🟡 仅 pTM 通过"
    elif ptm > 0.7:
        status = "🟠 接近 pTM 门槛"
    elif ptm > 0.5:
        status = "🔵 结构中等"
    else:
        status = "🔴 崩塌"
    
    rank += 1
    print(f"{rank:<5} {c['name'][:38]:<38} {c['scaffold'][:14]:<14} "
          f"{plddt:>6.1f} {chromo:>6.1f} {ptm:>6.4f} "
          f"{c['_gap_ptm']*100:>6.1f} {c['_gap_plddt']:>7.1f} {c['_gap_chromo']:>7.1f} {status:<20}")
    
    if rank >= 60:
        remaining = len(all_c) - rank
        print(f"  ... 还有 {remaining} 条 (全部低于 0.5 pTM)")
        break

print(f"\n{'='*80}")
print(f"统计:")
print(f"  通过全部结构门槛: {grades['PASS_ALL']} 条")
print(f"  通过 pTM>0.75: {grades['PASS_pTM']} 条")
print(f"  共评估: {len(all_c)} 条")
print(f"{'='*80}")

print(f"\n{'='*80}")
print(f"差距分析 (最接近的候选):")
print(f"{'='*80}")
for c in all_c[:5]:
    ptm = c.get("ptm", 0) or 0
    plddt = c.get("plddt_mean", 0) or 0
    chromo = c.get("plddt_chromo_region", plddt) or plddt
    print(f"\n  {c['name']} ({c['scaffold']})")
    print(f"    pTM:      {ptm:.4f} (需要 > 0.75, 差 {max(0,0.75-ptm)*100:.1f}%)")
    print(f"    pLDDT:    {plddt:.1f} (需要 > 80.0, 差 {max(0,80-plddt):.1f})")
    print(f"    Chromo:   {chromo:.1f} (需要 > 85.0, 差 {max(0,85-chromo):.1f})")
    print(f"    突变数:   {c['n_muts']}")
    print(f"    门槛满足: {'pTM' if ptm>0.75 else ''} {'pLDDT' if plddt>80 else ''} {'Chromo' if chromo>85 else ''}")
