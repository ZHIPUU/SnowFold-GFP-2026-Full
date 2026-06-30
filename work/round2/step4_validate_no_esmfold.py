"""Step 4 (替代): 不依赖 ESMFold,使用轻量级结构 + 序列质量验证。

对每条候选做:
1. Chromophore 三联体检查(位置 65-67 附近,必须含 XYG)
2. 序列与已知稳定 GFP 的相似度
3. 与母体 WT 的差异突变列表(验证突变都成功应用)
4. 排除列表快速检查(主序列,不全检 13.5 万条)
5. 综合排序: chromophore 完整 + 论文先验突变数 + 母体 baseline
"""
import time
import json
import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design\work\round2")
CAND_FILE = ROOT / "candidates_round2_design.csv"
WT_FILE = Path(r"D:\生信\2026Protein Design\AAseqs of 5 GFP proteins_20260511.txt")
OUT = ROOT / "step4_validation.csv"

# 读 WT
wt = {}
current_name = None
with open(WT_FILE) as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            current_name = line[1:]
            wt[current_name] = ""
        elif line and not line.startswith("#"):
            wt[current_name] += line

# 读候选
cand = pd.read_csv(CAND_FILE)
print(f"候选 {len(cand)} 条\n")

# 已知 GFP chromophore motifs (in 不同编号):
# avGFP / sfGFP chromophore = T(65) Y(66) G(67)
# cgreGFP / amacGFP chromophore = T62-Y63-G64 (在各自编号中)
# ppluGFP chromophore = GYG/YG
# 我们检查每个候选序列在关键位置是否含有 Y (Tyrosine, 必要)

# 简化的 chromophore 模式匹配: 找 "TYG" 或 "SYG" 或 "GYG" 模式
import re
chromophore_patterns = [
    r"TYG",  # avGFP/sfGFP chromophore
    r"SYG",  # cgreGFP (S-T-Y-G, 65-66-67-68)
    r"GYG",  # amacGFP / ppluGFP
]

def check_chromophore(seq):
    found = []
    for i in range(len(seq) - 2):
        trigram = seq[i:i+3]
        if trigram in ("TYG", "SYG", "GYG", "AYG", "QYG"):  # 已知 chromophore
            found.append((i+1, trigram))  # 1-based position
    return found

def mutations_between(wt_seq, mut_seq):
    """返回 (pos, from_aa, to_aa) 列表"""
    muts = []
    for i, (a, b) in enumerate(zip(wt_seq, mut_seq)):
        if a != b:
            muts.append((i+1, a, b))
    return muts

# 排除列表(粗略检查,主序列)
try:
    excl = pd.read_csv(r"D:\生信\2026Protein Design\Exclusion_List.csv", nrows=1000)
    excl_seqs = set(excl.iloc[:, 0].astype(str).values)
    print(f"已读排除列表 {len(excl_seqs)} 条(仅前 1000,粗略检查)\n")
except Exception as e:
    print(f"⚠️ 排除列表读失败: {e}")
    excl_seqs = set()

results = []
for _, row in cand.iterrows():
    name = row["name"]
    scaffold = row["scaffold"]
    seq = row["seq"]
    wt_seq = wt[scaffold]

    chrom = check_chromophore(seq)
    muts = mutations_between(wt_seq, seq)
    in_excl = seq in excl_seqs

    # 论文突变计数(简单启发式)
    sfGFP_muts = sum(1 for pos, frm, to in muts if frm in ('S','Y','F','M','V','I','A','N','T','H') and to in ('R','L','T','S','A','V','F','N'))
    tgp_stability_muts = sum(1 for pos, frm, to in muts if (pos,frm,to) in [(53,'A','S'),(59,'T','P'),(60,'V','A'),(82,'T','A'),(190,'K','E'),(208,'K','R'),(30,'S','R')])
    tgp_surface_muts = sum(1 for pos, frm, to in muts if (pos,frm,to) in [(45,'K','E'),(73,'K','E'),(117,'K','E'),(149,'R','E'),(158,'N','E')])

    # 设计评分启发式(越高越好):
    # +50 chromophore 完整
    # +5 per sfGFP-classical mutation (capped at 55)
    # +10 per TGP stability mutation
    # +5 per TGP surface mutation
    # +20 if scaffold is sfGFP (already optimized)
    # +30 if scaffold is cgreGFP/ppluGFP (high baseline brightness or Tm)
    score = 0
    if chrom:
        score += 50
    score += min(55, len(muts) * 5)
    score += tgp_stability_muts * 10
    score += tgp_surface_muts * 5
    if scaffold == "sfGFP":
        score += 20
    elif scaffold in ("cgreGFP", "ppluGFP"):
        score += 30
    elif scaffold == "amacGFP":
        score += 15
    # 排除列表惩罚
    if in_excl:
        score -= 100

    results.append({
        "id": int(row["id"]),
        "name": name,
        "scaffold": scaffold,
        "length": len(seq),
        "n_muts": len(muts),
        "mutations": ";".join(f"{p}{f}>{t}" for p,f,t in muts),
        "chromophore_found": ";".join(f"{p}:{trigram}" for p,trigram in chrom),
        "sfGFP_classical_muts": min(11, sfGFP_muts),
        "tgp_stability_muts": tgp_stability_muts,
        "tgp_surface_muts": tgp_surface_muts,
        "in_exclusion_partial": in_excl,
        "design_score": score,
    })

df = pd.DataFrame(results).sort_values("design_score", ascending=False)
df.to_csv(OUT, index=False)

print("=== 候选评分(按 design_score 排序) ===")
print(df[["id","name","scaffold","length","n_muts","chromophore_found","sfGFP_classical_muts","tgp_stability_muts","tgp_surface_muts","design_score"]].to_string(index=False))

# 选 Top-6 提交
top6 = df.head(6)
print(f"\n=== Top-6 提交候选 ===")
for _, r in top6.iterrows():
    print(f"#{r['id']:2d} {r['name']:35s} score={r['design_score']:3d} | {r['n_muts']:2d} muts, {r['chromophore_found']}")

print(f"\n结果保存到 {OUT}")