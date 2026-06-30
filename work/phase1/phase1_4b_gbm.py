"""
Phase 1.4b: 仅训练 LightGBM/sklearn GBM, 不做组合搜索
"""
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import r2_score
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
print(f"Total: {len(df)}")

# 训练 GBM
results = {}
for t in ["avGFP", "amacGFP", "cgreGFP", "ppluGFP"]:
    print(f"\n=== {t} ===")
    sub = df[(df["GFP type"] == t) & (df["n_mut"] <= 6)].reset_index(drop=True)
    n = len(sub)
    print(f"  n={n}")
    vocab = set()
    for pm in sub["parsed_muts"]:
        for k in pm:
            vocab.add(k)
    vocab = sorted(vocab)
    vocab_idx = {v: i for i, v in enumerate(vocab)}
    v = len(vocab)
    print(f"  vocab={v}")

    X = np.zeros((n, v), dtype=np.float32)
    y = sub["Brightness"].values
    for i, pm in enumerate(sub["parsed_muts"]):
        for k in pm:
            if k in vocab_idx:
                X[i, vocab_idx[k]] = 1.0

    # 3-fold CV (节省时间)
    kf = KFold(n_splits=3, shuffle=True, random_state=42)
    cv_r2 = []
    for fold, (tr, va) in enumerate(kf.split(X)):
        m = GradientBoostingRegressor(
            n_estimators=100, learning_rate=0.1, max_depth=5,
            min_samples_leaf=20, random_state=42
        )
        m.fit(X[tr], y[tr])
        pred = m.predict(X[va])
        r2 = r2_score(y[va], pred)
        cv_r2.append(r2)
        print(f"    fold {fold}: R2={r2:.3f}")
    print(f"  CV R2 = {np.mean(cv_r2):+.3f} (+/- {np.std(cv_r2):.3f})")
    results[t] = {"cv_r2_mean": float(np.mean(cv_r2)), "vocab_size": v, "n_train": n}

import json
with open(PHASE1 / "gbm_results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"\n[OK] Saved to {PHASE1 / 'gbm_results.json'}")
