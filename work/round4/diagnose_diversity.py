"""
Round 4 诊断: 检查 6 条候选的序列多样性
==================================
发现潜在问题: 多条候选可能高度相似, 同时失败的风险高
"""
import json
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
OUT = ROOT / "work" / "round4"

with open(OUT / "final_6_round4.json", encoding="utf-8") as f:
    final_6 = json.load(f)

# 两两对比汉明距离
print("=" * 80)
print("Round 4 提交 6 条候选 — 两两汉明距离矩阵:")
print("=" * 80)
print(f"{'':<6}", end="")
for i in range(6):
    print(f"{i+1:>5}", end="")
print()

for i, c1 in enumerate(final_6):
    print(f"Seq{i+1:<3}", end="  ")
    for j, c2 in enumerate(final_6):
        if i == j:
            print(f"{'-':>5}", end="")
        else:
            # 长度可能不同, 截断
            min_len = min(len(c1["seq"]), len(c2["seq"]))
            d = sum(1 for a, b in zip(c1["seq"][:min_len], c2["seq"][:min_len]) if a != b)
            if abs(len(c1["seq"]) - len(c2["seq"])) > 0:
                d += abs(len(c1["seq"]) - len(c2["seq"]))
            print(f"{d:>5}", end="")
    print(f"  ← {c1['name']}")

# 名称对照
print("\n候选对应:")
for i, c in enumerate(final_6, 1):
    print(f"  Seq {i}: {c['name']:<35} ({c['scaffold']}, {c['n_muts']} mut, Tm={c['expected_tm']})")

# 同骨架统计
from collections import Counter
print("\n骨架分布:")
for s, n in Counter(c["scaffold"] for c in final_6).most_common():
    print(f"  {s}: {n} 条")

# 多样性评分
print("\n问题诊断:")
sf_count = Counter(c["scaffold"] for c in final_6)
if sf_count.get("sfGFP", 0) >= 4:
    print(f"  ⚠️ sfGFP 占了 {sf_count['sfGFP']}/6 — 系统性风险高 (一坏全坏)")

# 序列相似度统计 (汉明距离 < 5 视为高相似)
n_similar = 0
for i in range(6):
    for j in range(i+1, 6):
        c1, c2 = final_6[i], final_6[j]
        if c1["scaffold"] == c2["scaffold"]:
            min_len = min(len(c1["seq"]), len(c2["seq"]))
            d = sum(1 for a, b in zip(c1["seq"][:min_len], c2["seq"][:min_len]) if a != b)
            if d < 5:
                n_similar += 1
                print(f"  ⚠️ Seq {i+1} vs Seq {j+1}: 汉明距离 {d} (过度相似)")

# Tm 分布
print("\nTm 分布:")
tm_counts = Counter(c["expected_tm"] for c in final_6)
for tm, n in sorted(tm_counts.items()):
    print(f"  Tm {tm}°C: {n} 条")

# 突变数分布
print("\n突变数分布:")
mut_counts = Counter(c["n_muts"] for c in final_6)
for m, n in sorted(mut_counts.items()):
    print(f"  {m} mut: {n} 条")
