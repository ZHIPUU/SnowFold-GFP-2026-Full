"""
Phase 1.5: 高效组合搜索 (纯加性模型, 向量化)
================================================
跳过 GBM, 直接用 Ridge 加性模型做组合搜索
- 限制: 2-6 突变
- 向量化计算, 避免 Python loop 卡顿
"""
import pickle
import sys
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

def parse_mut(s):
    return []
sys.modules['__main__'].parse_mut = parse_mut

WORK = Path(r"D:\生信\2026Protein Design\work")
PHASE1 = WORK / "phase1"

with open(PHASE1 / "additive_models_v2.pkl", "rb") as f:
    models = pickle.load(f)

with open(PHASE1 / "phase1_cache_v2.pkl", "rb") as f:
    cache = pickle.load(f)

wt_seqs = cache["wt_seqs"]
wt_brightness = cache["wt_brightness"]

# Exclusion List
print("Loading exclusion list...")
excl_df = pd.read_csv(r"D:\生信\2026Protein Design\Exclusion_List.csv")
excl_seqs = set(excl_df["Sequence"].values)
print(f"  {len(excl_seqs)} sequences in exclusion list")

# ---------- 组合搜索 ----------
def apply_mutations(wt_seq, muts):
    """应用突变列表到 WT 序列"""
    seq = list(wt_seq)
    for pos, wt_aa, new_aa in muts:
        idx = pos - 1
        if idx < 0 or idx >= len(seq):
            return None
        seq[idx] = new_aa
    return "".join(seq)

def check_constraints(seq, target_type):
    """检查长度、起始字符等约束"""
    if not seq or len(seq) < 220 or len(seq) > 250:
        return False
    if not seq.startswith("M"):
        return False
    if "*" in seq:
        return False
    # 排除列表
    if seq in excl_seqs:
        return False
    return True

# 收集所有候选 (按 type 分组)
top_per_type = {}
for t in ["avGFP", "amacGFP", "cgreGFP", "ppluGFP"]:
    print(f"\n=== {t} ===")
    m = models[t]
    intercept = m["intercept"]
    effects = m["effects"]

    # 排序 positive 效应
    pos_eff = sorted(
        [(k, v) for k, v in effects.items() if v > 0],
        key=lambda x: -x[1]
    )
    # top-30 阳性突变(经验: 多了搜索空间爆炸, 少了覆盖度不够)
    top_pos = pos_eff[:30]
    print(f"  Top {len(top_pos)} positive mutations")
    if len(top_pos) < 5:
        continue

    keys = []
    eff_dict = {}
    for k, v in top_pos:
        # k 是 (pos, wt_aa, new_aa) 或 (pos, new_aa)
        if len(k) == 3:
            pos, wt_aa, new_aa = k
        else:
            pos, new_aa = k
        keys.append((pos, new_aa))
        eff_dict[(pos, new_aa)] = v
    wt_seq = wt_seqs[t]

    candidates = []

    # 2-突变组合: C(30,2) = 435
    print(f"  Searching 2-mut...")
    for combo in combinations(keys, 2):
        if combo[0][0] == combo[1][0]:
            continue
        muts = []
        for pos, new_aa in combo:
            wt_aa = wt_seq[pos - 1] if pos - 1 < len(wt_seq) else "?"
            muts.append((pos, wt_aa, new_aa))
        seq = apply_mutations(wt_seq, muts)
        if not check_constraints(seq, t):
            continue
        pred = intercept + sum(eff_dict[k] for k in combo)
        candidates.append({"type": t, "n_mut": 2, "mutations": list(combo), "seq": seq, "pred": pred})

    print(f"    {len(candidates)} valid 2-mut candidates")

    # 3-突变组合: C(30,3) = 4060
    print(f"  Searching 3-mut...")
    cnt_3 = 0
    for combo in combinations(keys, 3):
        pos_list = [k[0] for k in combo]
        if len(set(pos_list)) != 3:
            continue
        cnt_3 += 1
        muts = []
        for pos, new_aa in combo:
            wt_aa = wt_seq[pos - 1] if pos - 1 < len(wt_seq) else "?"
            muts.append((pos, wt_aa, new_aa))
        seq = apply_mutations(wt_seq, muts)
        if not check_constraints(seq, t):
            continue
        pred = intercept + sum(eff_dict[k] for k in combo)
        candidates.append({"type": t, "n_mut": 3, "mutations": list(combo), "seq": seq, "pred": pred})
    print(f"    searched {cnt_3} 3-mut combos, {len([c for c in candidates if c['n_mut']==3])} valid")

    # 4-突变组合: C(30,4) = 27405
    print(f"  Searching 4-mut...")
    cnt_4 = 0
    for combo in combinations(keys, 4):
        pos_list = [k[0] for k in combo]
        if len(set(pos_list)) != 4:
            continue
        cnt_4 += 1
        muts = []
        for pos, new_aa in combo:
            wt_aa = wt_seq[pos - 1] if pos - 1 < len(wt_seq) else "?"
            muts.append((pos, wt_aa, new_aa))
        seq = apply_mutations(wt_seq, muts)
        if not check_constraints(seq, t):
            continue
        pred = intercept + sum(eff_dict[k] for k in combo)
        candidates.append({"type": t, "n_mut": 4, "mutations": list(combo), "seq": seq, "pred": pred})
    print(f"    searched {cnt_4} 4-mut combos, {len([c for c in candidates if c['n_mut']==4])} valid")

    # 5-突变组合: 贪心 (从 top 4-mut 加一个最优 1-mut)
    print(f"  Searching 5-mut (greedy from top 4-mut)...")
    four_mut_pool = [c for c in candidates if c["n_mut"] == 4]
    four_mut_pool.sort(key=lambda x: -x["pred"])
    cnt_5 = 0
    for base in four_mut_pool[:50]:  # 取 top 50 4-mut
        base_pos = set(k[0] for k in base["mutations"])
        for k in keys:
            if k[0] in base_pos:
                continue
            new_combo = base["mutations"] + [k]
            muts = []
            for pos, new_aa in new_combo:
                wt_aa = wt_seq[pos - 1] if pos - 1 < len(wt_seq) else "?"
                muts.append((pos, wt_aa, new_aa))
            seq = apply_mutations(wt_seq, muts)
            if not check_constraints(seq, t):
                continue
            pred = intercept + sum(eff_dict[kk] for kk in new_combo)
            cnt_5 += 1
            candidates.append({"type": t, "n_mut": 5, "mutations": list(new_combo), "seq": seq, "pred": pred})
    print(f"    {cnt_5} 5-mut candidates")

    # 6-突变: 贪心从 top 5-mut
    print(f"  Searching 6-mut (greedy from top 5-mut)...")
    five_mut_pool = [c for c in candidates if c["n_mut"] == 5]
    five_mut_pool.sort(key=lambda x: -x["pred"])
    cnt_6 = 0
    for base in five_mut_pool[:30]:
        base_pos = set(k[0] for k in base["mutations"])
        for k in keys:
            if k[0] in base_pos:
                continue
            new_combo = base["mutations"] + [k]
            muts = []
            for pos, new_aa in new_combo:
                wt_aa = wt_seq[pos - 1] if pos - 1 < len(wt_seq) else "?"
                muts.append((pos, wt_aa, new_aa))
            seq = apply_mutations(wt_seq, muts)
            if not check_constraints(seq, t):
                continue
            pred = intercept + sum(eff_dict[kk] for kk in new_combo)
            cnt_6 += 1
            candidates.append({"type": t, "n_mut": 6, "mutations": list(new_combo), "seq": seq, "pred": pred})
    print(f"    {cnt_6} 6-mut candidates")

    # 排序 + 去重 (按 seq)
    seen = set()
    unique = []
    for c in candidates:
        if c["seq"] in seen:
            continue
        seen.add(c["seq"])
        unique.append(c)
    candidates = sorted(unique, key=lambda x: -x["pred"])
    top_per_type[t] = candidates
    print(f"  Total unique valid: {len(candidates)}")
    print(f"  Top 10 by predicted brightness:")
    for c in candidates[:10]:
        mut_str = ":".join(f"{wt_seq[k[0]-1]}{k[0]}{k[1]}" for k in c["mutations"])
        wt_pred = wt_brightness[t]
        delta = c["pred"] - wt_pred
        rel = 10 ** delta
        print(f"    {c['n_mut']}-mut: pred={c['pred']:.3f} (Δ={delta:+.3f}, ~{rel:.2f}× WT)  {mut_str}")

# ---------- 全局汇总 ----------
all_candidates = []
for t, cands in top_per_type.items():
    all_candidates.extend(cands)
print(f"\n=== Total candidates: {len(all_candidates)} ===")

# 按 predicted brightness 排序, 跨类型
all_candidates.sort(key=lambda x: -x["pred"])
print("Top 30 overall:")
for c in all_candidates[:30]:
    mut_str = ":".join(f"{wt_seqs[c['type']][k[0]-1]}{k[0]}{k[1]}" for k in c["mutations"])
    wt_pred = wt_brightness[c["type"]]
    delta = c["pred"] - wt_pred
    rel = 10 ** delta
    print(f"  {c['type']:8s} {c['n_mut']}-mut: pred={c['pred']:.3f} (Δ={delta:+.3f}, ~{rel:.2f}× WT)  {mut_str}")

# 保存
df_out = pd.DataFrame([{
    "type": c["type"],
    "n_mut": c["n_mut"],
    "mutations": str(c["mutations"]),
    "mut_str": ":".join(f"{wt_seqs[c['type']][k[0]-1]}{k[0]}{k[1]}" for k in c["mutations"]),
    "pred_brightness": c["pred"],
    "pred_relative": 10 ** (c["pred"] - wt_brightness[c["type"]]),
    "seq": c["seq"],
} for c in all_candidates])
df_out.to_csv(PHASE1 / "top_candidates.csv", index=False)
print(f"\n[OK] Saved {len(df_out)} candidates to {PHASE1 / 'top_candidates.csv'}")
print("\n=== Phase 1.5 DONE ===")