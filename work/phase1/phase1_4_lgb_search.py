"""
Phase 1.4: 组合搜索 + LightGBM 对比
================================================
目标:
  1) 训练 LightGBM 模型 (非线性) 与 Ridge (加性) 对比
  2) 用加性模型进行"组合搜索" - 枚举 top-N 高产组合
  3) 产出 Top-100 候选供后续 Phase 2/3 评估
  4) 排除列表预筛 (Exclusion_List)
"""
import pickle
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import KFold

def parse_mut(s):
    return []
sys.modules['__main__'].parse_mut = parse_mut

WORK = Path(r"D:\生信\2026Protein Design\work")
PHASE1 = WORK / "phase1"

with open(PHASE1 / "phase1_cache_v2.pkl", "rb") as f:
    cache = pickle.load(f)

wt_seqs = cache["wt_seqs"]
wt_brightness = cache["wt_brightness"]
df = cache["df"].copy()

AA = "ACDEFGHIKLMNPQRSTVWY"
def parse_mut_v2(s):
    if s == "WT":
        return []
    out = []
    for tok in s.split(":"):
        i = 0
        while i < len(tok) and tok[i] in AA:
            i += 1
        j = len(tok) - 1
        while j >= 0 and tok[j] in AA:
            j -= 1
        if i == 0 or j == len(tok) - 1 or i > j:
            continue
        wt_aa = tok[:i]
        try:
            pos = int(tok[i:j + 1])
        except ValueError:
            continue
        new_aa = tok[j + 1:]
        if len(wt_aa) != 1 or len(new_aa) != 1:
            continue
        if wt_aa not in AA or new_aa not in AA:
            continue
        out.append((pos, wt_aa, new_aa))
    return out

df["parsed_muts"] = df["aaMutations"].apply(parse_mut_v2)
df["n_mut"] = df["parsed_muts"].apply(len)
df = df[df["n_mut"] > 0].reset_index(drop=True)
print(f"Total mutated rows: {len(df)}")

# 加载已有加性模型
with open(PHASE1 / "additive_models_v2.pkl", "rb") as f:
    additive_models = pickle.load(f)

# ---------- 1. 训练 LightGBM 模型 ----------
print("\n=== Training LightGBM models ===")
try:
    import lightgbm as lgb
    HAS_LGB = True
except ImportError:
    print("  LightGBM not installed, using sklearn GBM instead")
    from sklearn.ensemble import GradientBoostingRegressor
    HAS_LGB = False

lgb_models = {}
for t in ["avGFP", "amacGFP", "cgreGFP", "ppluGFP"]:
    print(f"\n  {t}:")
    sub = df[df["GFP type"] == t].reset_index(drop=True)
    sub = sub[sub["n_mut"] <= 8].reset_index(drop=True)
    n = len(sub)
    print(f"    n={n}")

    # 用全部 (pos, new_aa) 作为 vocab
    vocab = set()
    for pm in sub["parsed_muts"]:
        for k in pm:
            vocab.add(k)
    vocab = sorted(vocab)
    vocab_idx = {v: i for i, v in enumerate(vocab)}
    v = len(vocab)
    print(f"    vocab={v}")

    X = np.zeros((n, v), dtype=np.float32)
    y = sub["Brightness"].values
    for i, pm in enumerate(sub["parsed_muts"]):
        for k in pm:
            if k in vocab_idx:
                X[i, vocab_idx[k]] = 1.0

    # 5-fold CV
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    cv_r2 = []
    if HAS_LGB:
        for fold, (tr, va) in enumerate(kf.split(X)):
            m = lgb.LGBMRegressor(
                n_estimators=300, learning_rate=0.05, num_leaves=63,
                min_child_samples=20, reg_alpha=0.1, reg_lambda=0.1,
                n_jobs=-1, verbose=-1
            )
            m.fit(X[tr], y[tr])
            pred = m.predict(X[va])
            r2 = r2_score(y[va], pred)
            cv_r2.append(r2)
            print(f"      fold {fold}: R2={r2:.3f}")
        print(f"    CV R2 = {np.mean(cv_r2):+.3f} (+/- {np.std(cv_r2):.3f})")
        # 全数据
        m_full = lgb.LGBMRegressor(
            n_estimators=300, learning_rate=0.05, num_leaves=63,
            min_child_samples=20, reg_alpha=0.1, reg_lambda=0.1,
            n_jobs=-1, verbose=-1
        )
        m_full.fit(X, y)
        lgb_models[t] = {"model": m_full, "vocab": vocab, "vocab_idx": vocab_idx, "cv_r2": float(np.mean(cv_r2))}
    else:
        # sklearn GBM
        for fold, (tr, va) in enumerate(kf.split(X)):
            m = GradientBoostingRegressor(n_estimators=200, learning_rate=0.05, max_depth=6, random_state=42)
            m.fit(X[tr], y[tr])
            pred = m.predict(X[va])
            r2 = r2_score(y[va], pred)
            cv_r2.append(r2)
            print(f"      fold {fold}: R2={r2:.3f}")
        print(f"    CV R2 = {np.mean(cv_r2):+.3f} (+/- {np.std(cv_r2):.3f})")
        m_full = GradientBoostingRegressor(n_estimators=200, learning_rate=0.05, max_depth=6, random_state=42)
        m_full.fit(X, y)
        lgb_models[t] = {"model": m_full, "vocab": vocab, "vocab_idx": vocab_idx, "cv_r2": float(np.mean(cv_r2))}

# ---------- 2. 加载 Exclusion List ----------
print("\n=== Loading Exclusion List ===")
exclusion_df = pd.read_csv(r"D:\生信\2026Protein Design\Exclusion_List.csv")
exclusion_seqs = set(exclusion_df["Sequence"].values)
print(f"  Loaded {len(exclusion_seqs)} exclusion sequences")

# ---------- 3. 组合搜索 ----------
print("\n=== Combinatorial Search ===")
top_candidates = []
N_TOP_PER_TYPE = 30  # 每个 GFP type 选 top-30 组合
N_MUT_RANGE = [2, 3, 4, 5, 6, 7]  # 突变数范围

for t in ["avGFP", "amacGFP", "cgreGFP", "ppluGFP"]:
    print(f"\n--- {t} ---")
    m = additive_models[t]
    intercept = m["intercept"]
    effects = m["effects"]

    # 排序所有 positive 效应
    pos_eff = [(k, v) for k, v in effects.items() if v > 0]
    pos_eff.sort(key=lambda x: -x[1])
    # 取 top-50 阳性突变
    top_pos = pos_eff[:50]
    print(f"  Top {len(top_pos)} positive mutations considered")
    if len(top_pos) < 5:
        continue

    # 枚举所有 2-7 突变组合 (限制)
    # 2-突变: C(50,2)=1225
    # 3-突变: C(50,3)=19600
    # 4-突变: C(50,4)=230300
    # 5-突变: 2118760 (太多了)
    # 我们只枚举到 4-突变,5+ 用贪心
    candidates = []
    keys = [k for k, _ in top_pos]
    effects_dict = {k: v for k, v in top_pos}

    for n_mut in N_MUT_RANGE:
        if n_mut <= 4:
            for combo in combinations(keys, n_mut):
                # 检查位置不重复
                positions = [k[0] for k in combo]
                if len(set(positions)) != n_mut:
                    continue
                pred = intercept + sum(effects_dict[k] for k in combo)
                candidates.append((combo, pred, n_mut))
        else:
            # 贪心: 起始于 top 组合,逐步加最优突变
            # 简单实现: 从所有 4-突变组合中,尝试替换/添加 1 个
            base_combos = [c for c in candidates if c[2] == 4]
            base_combos.sort(key=lambda x: -x[1])
            for base, _, _ in base_combos[:50]:  # top 50 of 4-mut
                base_set = set(base)
                for k in keys:
                    if k in base_set:
                        continue
                    if k[0] in [kk[0] for kk in base_set]:
                        continue
                    new_combo = base + (k,)
                    pred = intercept + sum(effects_dict[kk] for kk in new_combo)
                    candidates.append((new_combo, pred, n_mut))

    # 去重
    seen = set()
    unique = []
    for combo, pred, n in candidates:
        key = tuple(sorted(combo, key=lambda x: (x[0], x[1])))
        if key in seen:
            continue
        seen.add(key)
        unique.append((combo, pred, n))
    candidates = unique

    # 排序
    candidates.sort(key=lambda x: -x[1])
    print(f"  Total unique candidates: {len(candidates)}")
    print(f"  Top 5 by predicted brightness:")
    for c, p, n in candidates[:5]:
        mut_str = ":".join(f"{wt_seqs[t][k[0]-1]}{k[0]}{k[1]}" for k in c)
        print(f"    {n}-mut: pred={p:.3f}  {mut_str}")

    # 保留 top-N per type
    for c, p, n in candidates[:N_TOP_PER_TYPE]:
        top_candidates.append({
            "type": t,
            "n_mut": n,
            "mutations": list(c),
            "mut_str": ":".join(f"{wt_seqs[t][k[0]-1]}{k[0]}{k[1]}" for k in c),
            "pred_brightness": p,
        })

print(f"\n=== Total candidates: {len(top_candidates)} ===")
cand_df = pd.DataFrame(top_candidates)
cand_df.to_csv(PHASE1 / "top_candidates.csv", index=False)
print(f"Saved to {PHASE1 / 'top_candidates.csv'}")

# ---------- 4. 保存 LGB 模型 ----------
if lgb_models:
    with open(PHASE1 / "lgb_models.pkl", "wb") as f:
        # 不能 pickle LGB 完整对象(可能有兼容问题),只保存关键信息
        save_lgb = {t: {"vocab": m["vocab"], "cv_r2": m["cv_r2"]} for t, m in lgb_models.items()}
        pickle.dump(save_lgb, f)
    print(f"Saved LGB model info to {PHASE1 / 'lgb_models.pkl'}")

# ---------- 5. 报告 ----------
import json
report = {
    "models": {
        "additive_ridge": {t: {"cv_r2": m["cv_r2_mean"], "cv_rmse": m["cv_rmse_mean"]} for t, m in additive_models.items()},
    },
    "n_candidates_per_type": {t: sum(1 for c in top_candidates if c["type"] == t) for t in ["avGFP", "amacGFP", "cgreGFP", "ppluGFP"]},
}
if lgb_models:
    report["models"]["lightgbm"] = {t: {"cv_r2": m["cv_r2"]} for t, m in lgb_models.items()}
with open(PHASE1 / "phase1_4_report.json", "w") as f:
    json.dump(report, f, indent=2)
print(f"\n[OK] Phase 1.4 report saved")
print("\n=== Phase 1.4 DONE ===")
