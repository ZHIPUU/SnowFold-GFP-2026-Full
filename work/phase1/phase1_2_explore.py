"""
Phase 1.2: 14万条 GFP 突变-亮度数据深度分析
================================================
目标:
  1) 加载数据并整理
  2) 按 GFP 类型分组建模
  3) 计算单点突变的加性效应矩阵
  4) 识别关键 hotspot
  5) 评估突变数 vs 亮度的关系
  6) 估计 epistasis 信号
  7) 输出可读报告 + 缓存
"""
import json
import os
import pickle
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import openpyxl
import pandas as pd
import seaborn as sns
from tqdm import tqdm

WORK = Path(r"D:\生信\2026Protein Design\work")
PHASE1 = WORK / "phase1"
PHASE1.mkdir(parents=True, exist_ok=True)

XLSX = Path(r"D:\生信\2026Protein Design\GFP_data.xlsx")
WT_FILE = Path(r"D:\生信\2026Protein Design\AAseqs of 5 GFP proteins_20260511.txt")

# ---------- 1. 加载 WT 序列 ----------
wt_seqs = {}
with open(WT_FILE) as f:
    cur_name, cur_seq = None, []
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if cur_name:
                wt_seqs[cur_name] = "".join(cur_seq)
            cur_name = line[1:].split()[0]
            cur_seq = []
        elif not line.startswith("#"):
            cur_seq.append(line)
    if cur_name:
        wt_seqs[cur_name] = "".join(cur_seq)

print("WT sequences loaded:", {k: len(v) for k, v in wt_seqs.items()})

# ---------- 2. 加载 brightness 数据 ----------
wb = openpyxl.load_workbook(XLSX, read_only=True)
ws = wb["brightness"]
rows = list(ws.iter_rows(values_only=True))
header = rows[0]
df = pd.DataFrame(rows[1:], columns=header)
print(f"\nTotal rows: {len(df)}")
print(f"Columns: {list(df.columns)}")
print(f"df dtypes:\n{df.dtypes}")
print(f"df head:\n{df.head()}")

# ---------- 3. 基础分布 ----------
print("\n=== GFP type distribution ===")
print(df["GFP type"].value_counts())

print("\n=== Brightness by type (WT vs non-WT) ===")
for t in df["GFP type"].unique():
    sub = df[df["GFP type"] == t]
    wt_b = sub[sub["aaMutations"] == "WT"]["Brightness"].values
    nonwt = sub[sub["aaMutations"] != "WT"]["Brightness"].values
    print(f"  {t}: WT={wt_b}, n_nonWT={len(nonwt)}, max={nonwt.max():.3f}, "
          f"mean={nonwt.mean():.3f}, top5={sorted(nonwt, reverse=True)[:5]}")

# ---------- 4. 突变数 vs 亮度 ----------
print("\n=== Mutation count vs brightness ===")
df["n_mut"] = df["aaMutations"].apply(lambda x: 0 if x == "WT" else len(x.split(":")))
for t in df["GFP type"].unique():
    sub = df[df["GFP type"] == t]
    print(f"\n  {t}:")
    for n in sorted(sub["n_mut"].unique())[:10]:
        sb = sub[sub["n_mut"] == n]["Brightness"]
        if len(sb) > 0:
            print(f"    n_mut={n}: count={len(sb)}, mean={sb.mean():.3f}, "
                  f"max={sb.max():.3f}, std={sb.std():.3f}")

# ---------- 5. 解析单点突变 ----------
# 格式: <WT_AA><pos><new_AA>  例如 A109D
print("\n=== Parsing single-point mutations ===")
AA = "ACDEFGHIKLMNPQRSTVWY"

def parse_mut(s):
    if s == "WT":
        return []
    out = []
    for tok in s.split(":"):
        # 找到第一个数字位置
        i = 0
        while i < len(tok) and tok[i] in AA:
            i += 1
        j = len(tok) - 1
        while j >= 0 and tok[j] in AA:
            j -= 1
        if i == 0 or j == len(tok) - 1 or i > j:
            continue
        wt_aa = tok[:i]
        pos = int(tok[i:j + 1])
        new_aa = tok[j + 1:]
        if len(wt_aa) != 1 or len(new_aa) != 1:
            continue
        if wt_aa not in AA or new_aa not in AA:
            continue
        out.append((pos, wt_aa, new_aa))
    return out

# 抽取单点突变表
single_records = []
for _, r in df.iterrows():
    if r["aaMutations"] == "WT":
        continue
    muts = parse_mut(r["aaMutations"])
    if len(muts) == 1:
        pos, wt_aa, new_aa = muts[0]
        # 验证 WT aa 跟参考序列一致
        wt_seq = wt_seqs.get(r["GFP type"])
        if wt_seq and pos - 1 < len(wt_seq) and wt_seq[pos - 1] == wt_aa:
            single_records.append({
                "type": r["GFP type"],
                "pos": pos,
                "wt_aa": wt_aa,
                "new_aa": new_aa,
                "brightness": r["Brightness"],
            })
print(f"  Single-point records (validated against WT): {len(single_records)}")
single_df = pd.DataFrame(single_records)
print(f"  Single point records by type:\n{single_df['type'].value_counts()}")

# ---------- 6. 单点突变的 WT baseline (按 type) ----------
wt_brightness = {}
for t in df["GFP type"].unique():
    wt_brightness[t] = float(df[(df["GFP type"] == t) & (df["aaMutations"] == "WT")]["Brightness"].iloc[0])
print("\n=== WT brightness per type ===")
print(wt_brightness)

# 计算 delta brightness = log10(突变体 / WT)
single_df["delta_b"] = single_df.apply(
    lambda r: r["brightness"] - wt_brightness[r["type"]], axis=1
)
print(f"\n  delta_b stats:\n{single_df['delta_b'].describe()}")

# ---------- 7. 每个 type 的 hotspot (按位置聚合) ----------
print("\n=== Hotspot positions per type (mean delta_b for positions with >=3 obs) ===")
hotspots = {}
for t in single_df["type"].unique():
    sub = single_df[single_df["type"] == t]
    pos_stats = sub.groupby("pos").agg(
        n=("delta_b", "size"),
        mean_delta=("delta_b", "mean"),
        max_delta=("delta_b", "max"),
        min_delta=("delta_b", "min"),
    ).reset_index()
    pos_stats = pos_stats[pos_stats["n"] >= 3].sort_values("mean_delta", ascending=False)
    hotspots[t] = pos_stats
    print(f"\n  {t}: {len(pos_stats)} positions observed")
    print(f"    Top 10 by mean_delta_b:")
    for _, r in pos_stats.head(10).iterrows():
        wt_aa = wt_seqs[t][int(r["pos"]) - 1] if int(r["pos"]) - 1 < len(wt_seqs[t]) else "?"
        print(f"      pos {r['pos']} ({wt_aa}): n={r['n']}, "
              f"mean_d={r['mean_delta']:+.3f}, max_d={r['max_delta']:+.3f}, min_d={r['min_delta']:+.3f}")

# ---------- 8. 缓存 ----------
cache = {
    "wt_seqs": wt_seqs,
    "wt_brightness": wt_brightness,
    "single_df": single_df,
    "hotspots": hotspots,
    "df": df,
    "parse_mut": parse_mut,
}
with open(PHASE1 / "phase1_cache.pkl", "wb") as f:
    pickle.dump(cache, f)
print(f"\n[OK] Cache saved to {PHASE1 / 'phase1_cache.pkl'}")

# ---------- 9. 报告 ----------
report_lines = [
    "# Phase 1.2 数据深度分析报告",
    "",
    "## 1. 数据规模",
    f"- 总突变记录: {len(df)} 条",
    f"- 4 种 GFP 母体: avGFP / amacGFP / cgreGFP / ppluGFP",
    f"- 单点突变(已验证): {len(single_df)} 条",
    "",
    "## 2. WT 基线 (CFPS 初始亮度, log10 尺度)",
]
for t, b in wt_brightness.items():
    report_lines.append(f"- {t}: {b:.3f} (相对 = 1.0)")
report_lines += [
    "",
    "## 3. 突变数 vs 亮度 (以 avGFP 为例)",
]
for n in range(0, 8):
    sub = df[(df["GFP type"] == "avGFP") & (df["n_mut"] == n)]
    if len(sub):
        report_lines.append(
            f"- n_mut={n}: count={len(sub)}, mean={sub['Brightness'].mean():.3f}, max={sub['Brightness'].max():.3f}"
        )

report_lines += [
    "",
    "## 4. 各母体 Top-10 Hotspot (按 mean delta_brightness 排序, 至少 3 次观测)",
]
for t, hs in hotspots.items():
    report_lines.append(f"\n### {t}\n")
    for _, r in hs.head(10).iterrows():
        wt_aa = wt_seqs[t][int(r["pos"]) - 1] if int(r["pos"]) - 1 < len(wt_seqs[t]) else "?"
        report_lines.append(
            f"- pos {r['pos']:>3} ({wt_aa}): n={r['n']:>2}, "
            f"mean_d={r['mean_delta']:+.3f}, max_d={r['max_delta']:+.3f}"
        )

with open(PHASE1 / "phase1_report.md", "w", encoding="utf-8") as f:
    f.write("\n".join(report_lines))
print(f"[OK] Report saved to {PHASE1 / 'phase1_report.md'}")
print("\n=== Phase 1.2 DONE ===")
