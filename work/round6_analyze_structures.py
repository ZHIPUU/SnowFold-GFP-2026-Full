"""分析所有候选的结构质量"""
import json
from pathlib import Path

r5 = Path(r"D:\生信\2026Protein Design\work\round5")
with open(r5 / "evolvepro_scored.json", encoding="utf-8") as f:
    all_c = json.load(f)

# 按 pTM 排序
all_c.sort(key=lambda x: -(x.get("ptm",0) or 0))

print(f"{'#':<3} {'Name':<35} {'Scaffold':<16} {'pLDDT':>6} {'chromo':>6} {'pTM':>6} {'muts':>4}")
print("-" * 85)
for i, c in enumerate(all_c, 1):
    ptm = c.get("ptm",0) or 0
    plddt = c.get("plddt_mean",0) or 0
    chromo = c.get("plddt_chromo_region", plddt) or plddt
    print(f"{i:<3} {c['name'][:35]:<35} {c['scaffold'][:16]:<16} {plddt:>6.1f} {chromo:>6.1f} {ptm:>6.4f} {c['n_muts']:>4}")

# 统计
print(f"\n--- 统计 ---")
ptm_thresholds = [0.75, 0.7, 0.65, 0.6, 0.5]
for th in ptm_thresholds:
    n = sum(1 for c in all_c if (c.get("ptm",0) or 0) >= th)
    print(f"  pTM >= {th}: {n} 条")

chromo_thresholds = [80, 70, 60, 50]
for th in chromo_thresholds:
    n = sum(1 for c in all_c if (c.get("plddt_chromo_region", 0) or 0) >= th)
    print(f"  生色团 pLDDT >= {th}: {n} 条")

# 同时满足
for pt, ch in [(0.75, 80), (0.75, 70), (0.7, 70), (0.7, 60), (0.65, 60)]:
    n = sum(1 for c in all_c 
            if (c.get("ptm",0) or 0) >= pt 
            and (c.get("plddt_chromo_region", 0) or 0) >= ch)
    print(f"  pTM>={pt} + 生色团>={ch}: {n} 条")

# 列出 pTM 前10 的具体信息
print(f"\n--- pTM Top-10 详情 ---")
for i, c in enumerate(all_c[:10], 1):
    ptm = c.get("ptm",0) or 0
    plddt = c.get("plddt_mean",0) or 0
    chromo = c.get("plddt_chromo_region", plddt) or plddt
    muts = c["n_muts"]
    scaf = c["scaffold"]
    # 判断结构等级
    if ptm >= 0.75 and chromo >= 80:
        grade = "GOLD"
    elif ptm >= 0.7 and chromo >= 70:
        grade = "SILVER"
    elif ptm >= 0.65 and chromo >= 60:
        grade = "BRONZE"
    elif ptm >= 0.5:
        grade = "RISKY"
    else:
        grade = "COLLAPSED"
    print(f"  #{i:2d} {c['name'][:30]:<30} {scaf:<14} pLDDT={plddt:.1f} chromo={chromo:.1f} pTM={ptm:.4f} muts={muts:>3} [{grade}]")
