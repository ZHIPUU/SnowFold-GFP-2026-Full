"""
Phase 6: 按比赛规则评分预测
================================================
比赛规则:
  综合分 = (Finitial / Finitial_WT) × (Ffinal / Finitial)
  极低阈值: Finitial < 0.3 × Finitial_WT → 0 分

我们有:
  - Finitial/Finitial_WT: 加性 Ridge 模型预测 (R² 0.7+, 可信)
  - Ffinal/Finitial:     只能靠先验估算 (文献/已知超稳突变)
"""
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

WORK = Path(r"D:\生信\2026Protein Design\work")
PHASE1 = WORK / "phase1"
PHASE3 = WORK / "phase3"

# 加载
df = pd.read_csv(PHASE3 / "final_6_candidates.csv")
with open(PHASE1 / "phase1_cache_v2.pkl", "rb") as f:
    cache = pickle.load(f)
wt_seqs = cache["wt_seqs"]
wt_brightness = cache["wt_brightness"]
with open(PHASE1 / "additive_models_v2.pkl", "rb") as f:
    models = pickle.load(f)

# ---------- 文献 Tm 值 (用于热稳定性估算) ----------
# 各种 GFP 母体的 Tm (从文献,大致值)
TMs = {
    "avGFP": 62,    # 原始 GFP
    "sfGFP": 70,    # superfolder
    "amacGFP": 55,  # amac
    "cgreGFP": 65,  # cgre (copGFP 类似)
    "ppluGFP": 50,  # pplu
}

# 不同类型突变对 Tm 的影响 (基于 sfGFP 文献)
MUT_TM_IMPACT = {
    # sfGFP 经典超稳突变 (每加一个 +3-5°C)
    "S65T": 4,
    "S72A": 3,
    "V163A": 4,
    "T203I": 3,
    "S202D": 4,
    "A206V": 2,
    # cgreGFP hotspot
    "K167M": 5,
    # ppluGFP hotspot
    "L199H": 6,
    "T137S": 3,
    # amacGFP hotspot
    "I166T": 5,
}

def estimate_tm(type_name, mutations_str):
    """基于母体 Tm + 突变影响估算 Tm"""
    base_tm = TMs.get(type_name, 60)
    delta = 0
    for mut in mutations_str.split("/"):
        mut = mut.strip()
        for key, impact in MUT_TM_IMPACT.items():
            if mut.endswith(key) or mut == key:
                delta += impact
                break
    return base_tm + delta

def tm_to_finitial_ratio(tm, treatment_temp=72):
    """
    估算 72°C 热处理后剩余亮度比例 Ffinal/Finitial
    基于两态变性模型:
      fraction_folded = 1 / (1 + exp(-ΔH/R × (1/Tm - 1/T)))
    简化: 用 sigmoid 近似
    """
    if tm >= treatment_temp + 10:
        return 0.95
    elif tm >= treatment_temp:
        # 接近 Tm,残留 50%
        return 0.5 + 0.4 * (tm - treatment_temp) / 10
    elif tm >= treatment_temp - 10:
        # 低于 Tm,残留 10-50%
        return 0.1 + 0.4 * (tm - (treatment_temp - 10)) / 10
    else:
        # 远低于 Tm,基本变性
        return max(0.05, 0.1 * (tm / (treatment_temp - 10)))

# ---------- 计算评分 ----------
print("=" * 100)
print("按比赛规则评分预测")
print("=" * 100)
print(f"\n{'#':<3} {'Type':<10} {'nMut':<5} {'pred_b':<8} {'rel×':<6} "
      f"{'Tm_est':<7} {'Ffinal/Finitial':<18} {'综合分(乐观)':<12} {'综合分(中等)':<12} {'综合分(保守)':<12}")
print("-" * 130)

results = []
for _, row in df.iterrows():
    t = row["Type"]
    n_mut = row["n_mut"]
    muts = row["mutations"]

    # Finitial/Finitial_WT: 基于加性模型 (如果可计算)
    if pd.notna(row["pred_brightness"]):
        rel_b = 10 ** (row["pred_brightness"] - wt_brightness[t])
    else:
        # Seq 6 (sfGFP-classical): 加性模型不适用, 用文献近似
        # sfGFP 比 avGFP 亮度高 ~3-4×, 加 I152S 再 +0.528
        rel_b = 3.0 * 10 ** 0.528  # 估算 ~10× WT
        rel_b = min(rel_b, 10)  # cap at 10×

    # 30% 阈值检查
    if rel_b < 0.3:
        final_score = 0
        tm_est = None
        ratio_est = 0
    else:
        tm_est = estimate_tm(t, muts)
        # 三种情景
        ratio_optimistic = min(0.95, tm_to_finitial_ratio(tm_est + 5))
        ratio_moderate = tm_to_finitial_ratio(tm_est)
        ratio_conservative = max(0.05, tm_to_finitial_ratio(tm_est - 5))

        final_optimistic = rel_b * ratio_optimistic
        final_moderate = rel_b * ratio_moderate
        final_conservative = rel_b * ratio_conservative
    results.append({
        "Seq_ID": row["Seq_ID"],
        "Type": t,
        "n_mut": n_mut,
        "mutations": muts,
        "rel_b_initial": rel_b,
        "tm_est_C": tm_est,
        "ratio_opt": ratio_optimistic if rel_b >= 0.3 else 0,
        "ratio_mod": ratio_moderate if rel_b >= 0.3 else 0,
        "ratio_con": ratio_conservative if rel_b >= 0.3 else 0,
        "score_opt": final_optimistic if rel_b >= 0.3 else 0,
        "score_mod": final_moderate if rel_b >= 0.3 else 0,
        "score_con": final_conservative if rel_b >= 0.3 else 0,
    })

res_df = pd.DataFrame(results)
print(f"\n{'#':<3} {'Type':<10} {'nMut':<5} {'rel×':<6} "
      f"{'Tm°C':<6} {'F/F (opt)':<10} {'综合(乐观)':<10} {'综合(中等)':<10} {'综合(保守)':<10} {'Mutations'}")
print("-" * 130)
for _, r in res_df.iterrows():
    print(f"{int(r['Seq_ID']):<3} {r['Type']:<10} {int(r['n_mut']):<5} {r['rel_b_initial']:<6.2f} "
          f"{r['tm_est_C'] if r['tm_est_C'] else 'N/A':<6} {r['ratio_opt']:<10.2f} "
          f"{r['score_opt']:<10.2f} {r['score_mod']:<10.2f} {r['score_con']:<10.2f} {r['mutations'][:40]}")

# ---------- Top-1 排名 ----------
print("\n" + "=" * 80)
print("按比赛规则:Best Top-1 综合分作为打榜成绩 (中等情景)")
print("=" * 80)
top1_mod = res_df.sort_values("score_mod", ascending=False).iloc[0]
print(f"  Top-1 (中等): Seq {int(top1_mod['Seq_ID'])} ({top1_mod['Type']}, {top1_mod['n_mut']}-mut)")
print(f"    综合分 = {top1_mod['score_mod']:.2f}")
print(f"    Mutations: {top1_mod['mutations']}")
print(f"    Finitial/Finitial_WT ≈ {top1_mod['rel_b_initial']:.2f}×")
print(f"    Ffinal/Finitial ≈ {top1_mod['ratio_mod']:.2f}")
print(f"    Tm_est ≈ {top1_mod['tm_est_C']:.0f}°C")

top1_opt = res_df.sort_values("score_opt", ascending=False).iloc[0]
print(f"\n  Top-1 (乐观): Seq {int(top1_opt['Seq_ID'])} ({top1_opt['Type']}, {top1_opt['n_mut']}-mut)")
print(f"    综合分 = {top1_opt['score_opt']:.2f}")

top1_con = res_df.sort_values("score_con", ascending=False).iloc[0]
print(f"\n  Top-1 (保守): Seq {int(top1_con['Seq_ID'])} ({top1_con['Type']}, {top1_con['n_mut']}-mut)")
print(f"    综合分 = {top1_con['score_con']:.2f}")

# ---------- 最佳亮度奖 / 最佳热稳奖 ----------
print("\n" + "=" * 80)
print("独立单项奖")
print("=" * 80)
best_brightness = res_df.sort_values("rel_b_initial", ascending=False).iloc[0]
print(f"\n  最佳亮度 (相对 Finitial): Seq {int(best_brightness['Seq_ID'])} ({best_brightness['Type']})")
print(f"    Finitial/Finitial_WT ≈ {best_brightness['rel_b_initial']:.2f}×")
print(f"    Mutations: {best_brightness['mutations']}")

best_thermal = res_df.sort_values("ratio_mod", ascending=False).iloc[0]
print(f"\n  最佳热稳 (Ffinal/Finitial): Seq {int(best_thermal['Seq_ID'])} ({best_thermal['Type']})")
print(f"    Ffinal/Finitial ≈ {best_thermal['ratio_mod']:.2f}")
print(f"    Mutations: {best_thermal['mutations']}")
print(f"    Tm_est ≈ {best_thermal['tm_est_C']:.0f}°C")

# ---------- 与历史 30% 金奖阈值比较 ----------
print("\n" + "=" * 80)
print("与历史金奖阈值比较 (2024/2025 头部队伍水平)")
print("=" * 80)
print(f"  历史金奖预估阈值: 综合分 ≥ 1.0 进 30% 金奖")
print(f"  历史银奖预估阈值: 综合分 ≥ 0.5")
print(f"\n  本队 Top-1 综合分 (中等情景): {top1_mod['score_mod']:.2f}")
if top1_mod['score_mod'] >= 1.5:
    print(f"  >>> 强金奖级别 (>1.5)")
elif top1_mod['score_mod'] >= 1.0:
    print(f"  >>> 金奖级别 (1.0-1.5)")
elif top1_mod['score_mod'] >= 0.5:
    print(f"  >>> 银奖级别")
else:
    print(f"  >>> 需重新评估")

# 保存
res_df.to_csv(WORK / "phase6_scores.csv", index=False)
print(f"\n[OK] Saved to {WORK / 'phase6_scores.csv'}")