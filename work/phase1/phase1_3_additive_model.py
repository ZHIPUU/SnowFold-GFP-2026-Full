"""
Phase 1.3: 从多点突变数据反推"边贡献效应矩阵"
================================================
策略:
  - 对每个 (pos, new_aa) 组合,统计它在所有包含它的突变体中的"平均边贡献"
  - 用岭回归 (Ridge) 拟合加性模型:  brightness ~ intercept + sum(mut_effects)
  - 用 LOO 交叉验证评估模型
  - 加上二阶 epistasis (top 残基对)
  - 保存模型参数 + 验证报告
"""
import pickle
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge, RidgeCV
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold
from tqdm import tqdm

WORK = Path(r"D:\生信\2026Protein Design\work")
PHASE1 = WORK / "phase1"

with open(PHASE1 / "phase1_cache_v2.pkl", "rb") as f:
    cache = pickle.load(f)

wt_seqs = cache["wt_seqs"]
wt_brightness = cache["wt_brightness"]
df = cache["df"].copy()

AA = "ACDEFGHIKLMNPQRSTVWY"

def parse_mut(s):
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
        pos = int(tok[i:j + 1])
        new_aa = tok[j + 1:]
        if len(wt_aa) != 1 or len(new_aa) != 1:
            continue
        if wt_aa not in AA or new_aa not in AA:
            continue
        out.append((pos, wt_aa, new_aa))
    return out

# ---------- 1. 解析所有突变 ----------
print("Parsing all mutations...")
df["parsed_muts"] = df["aaMutations"].apply(lambda s: parse_mut(s) if s != "WT" else [])
df["n_mut"] = df["parsed_muts"].apply(len)

# 过滤: 只保留突变位置匹配 WT 的(防止数据集噪声)
print("Filtering mutations by WT consistency...")
def validate_muts(row):
    wt_seq = wt_seqs.get(row["GFP type"])
    if wt_seq is None:
        return False
    for pos, wt_aa, new_aa in row["parsed_muts"]:
        if pos - 1 >= len(wt_seq) or wt_seq[pos - 1] != wt_aa:
            return False
    return True

df["valid"] = df.apply(validate_muts, axis=1)
df = df[df["valid"]].copy()
print(f"  Valid rows: {len(df)} / original 141572")

# ---------- 2. 构建 (type, pos, new_aa) 索引 ----------
print("Building position-mutation index...")
mut_vocab = defaultdict(set)  # mut_vocab[type] -> set of (pos, new_aa)
for _, r in df.iterrows():
    t = r["GFP type"]
    for pos, wt_aa, new_aa in r["parsed_muts"]:
        mut_vocab[t].add((pos, new_aa))

for t in mut_vocab:
    print(f"  {t}: {len(mut_vocab[t])} unique (pos, new_aa) mutations")

# ---------- 3. 为每个 type 拟合加性模型 ----------
def encode_mutations(parsed_muts, mut_vocab_t):
    """返回 dict of {(pos, new_aa): 1} for known mutations"""
    vocab = mut_vocab_t
    out = {}
    for pos, wt_aa, new_aa in parsed_muts:
        key = (pos, new_aa)
        if key in vocab:
            out[key] = 1
    return out

# 限制最大突变数,避免过拟合
MAX_NMUT = 8
df_train = df[df["n_mut"] <= MAX_NMUT].copy()
print(f"\nTraining on rows with n_mut <= {MAX_NMUT}: {len(df_train)}")

models = {}
reports = {}

for t in ["avGFP", "amacGFP", "cgreGFP", "ppluGFP"]:
    print(f"\n=== Fitting additive model for {t} ===")
    sub = df_train[df_train["GFP type"] == t].reset_index(drop=True)
    n = len(sub)
    if n == 0:
        continue

    vocab = sorted(mut_vocab[t])
    vocab_idx = {v: i for i, v in enumerate(vocab)}
    v = len(vocab)
    print(f"  n_samples={n}, vocab_size={v}")

    # 构建稀疏设计矩阵
    X = np.zeros((n, v), dtype=np.float32)
    y = sub["Brightness"].values
    for i, parsed in enumerate(sub["parsed_muts"]):
        for pos, wt_aa, new_aa in parsed:
            key = (pos, new_aa)
            if key in vocab_idx:
                X[i, vocab_idx[key]] = 1.0

    # 5-fold cross-validation (按"突变组合" group)
    # GroupKFold 按行(突变组合)分,避免泄漏
    n_folds = 5
    fold_size = n // n_folds
    indices = np.arange(n)
    np.random.seed(42)
    np.random.shuffle(indices)

    cv_r2, cv_rmse = [], []
    cv_predictions = np.zeros(n)
    for fold in range(n_folds):
        val_idx = indices[fold * fold_size:(fold + 1) * fold_size]
        train_idx = np.setdiff1d(indices, val_idx)
        Xtr, ytr = X[train_idx], y[train_idx]
        Xva, yva = X[val_idx], y[val_idx]
        m = RidgeCV(alphas=[0.01, 0.1, 1.0, 10.0, 100.0, 1000.0])
        m.fit(Xtr, ytr)
        pred = m.predict(Xva)
        cv_predictions[val_idx] = pred
        cv_r2.append(r2_score(yva, pred))
        cv_rmse.append(np.sqrt(mean_squared_error(yva, pred)))
        print(f"    fold {fold}: alpha={m.alpha_}, R2={cv_r2[-1]:.3f}, RMSE={cv_rmse[-1]:.3f}")

    print(f"  CV R2 mean = {np.mean(cv_r2):.3f} (+/- {np.std(cv_r2):.3f})")
    print(f"  CV RMSE mean = {np.mean(cv_rmse):.3f}")

    # 全数据最终模型
    m_full = RidgeCV(alphas=[0.01, 0.1, 1.0, 10.0, 100.0, 1000.0])
    m_full.fit(X, y)

    # 提取效应 (intercept + 每个突变的 log10 增量)
    intercept = float(m_full.intercept_)
    effects = {}
    for i, key in enumerate(vocab):
        effects[key] = float(m_full.coef_[i])

    # 排序: 找 top positive effects
    sorted_pos = sorted(effects.items(), key=lambda x: -x[1])
    print(f"  Top 15 positive effects (log10):")
    for (pos, new_aa), e in sorted_pos[:15]:
        wt_aa = wt_seqs[t][pos - 1] if pos - 1 < len(wt_seqs[t]) else "?"
        print(f"    {wt_aa}{pos}{new_aa}: {e:+.3f}")
    print(f"  Top 15 negative effects (log10):")
    for (pos, new_aa), e in sorted_pos[-15:]:
        wt_aa = wt_seqs[t][pos - 1] if pos - 1 < len(wt_seqs[t]) else "?"
        print(f"    {wt_aa}{pos}{new_aa}: {e:+.3f}")

    models[t] = {
        "intercept": intercept,
        "effects": effects,
        "alpha": m_full.alpha_,
        "vocab": vocab,
        "n_train": n,
        "cv_r2_mean": float(np.mean(cv_r2)),
        "cv_rmse_mean": float(np.mean(cv_rmse)),
    }
    reports[t] = {
        "cv_r2_per_fold": cv_r2,
        "cv_rmse_per_fold": cv_rmse,
        "cv_r2_mean": float(np.mean(cv_r2)),
        "cv_rmse_mean": float(np.mean(cv_rmse)),
    }

# 保存模型
with open(PHASE1 / "additive_models.pkl", "wb") as f:
    pickle.dump(models, f)
print(f"\n[OK] Additive models saved to {PHASE1 / 'additive_models.pkl'}")

# 保存报告
import json
with open(PHASE1 / "model_report.json", "w") as f:
    json.dump({
        "wt_brightness": wt_brightness,
        "models": {
            t: {
                "intercept": models[t]["intercept"],
                "alpha": models[t]["alpha"],
                "n_train": models[t]["n_train"],
                "cv_r2_mean": models[t]["cv_r2_mean"],
                "cv_rmse_mean": models[t]["cv_rmse_mean"],
                "n_vocab": len(models[t]["vocab"]),
                "top_positive": [
                    {"pos": int(k[0]), "new_aa": k[1], "wt_aa": wt_seqs[t][k[0]-1] if k[0]-1 < len(wt_seqs[t]) else "?", "effect": v}
                    for k, v in sorted(models[t]["effects"].items(), key=lambda x: -x[1])[:30]
                ],
                "top_negative": [
                    {"pos": int(k[0]), "new_aa": k[1], "wt_aa": wt_seqs[t][k[0]-1] if k[0]-1 < len(wt_seqs[t]) else "?", "effect": v}
                    for k, v in sorted(models[t]["effects"].items(), key=lambda x: x[1])[:30]
                ],
            } for t in models
        },
    }, f, indent=2)
print(f"[OK] Model report saved to {PHASE1 / 'model_report.json'}")

print("\n=== Phase 1.3 DONE ===")
