"""
Phase 1.3 v2: 放宽数据过滤,使用全部数据
================================================
策略:
  - 不再强制要求"突变描述中的 WT AA 跟参考序列完全一致"
  - 因为数据集可能基于"中间变体"做的扫描
  - 用全部 14 万条数据,但把单点 / 多点突变都纳入
  - 用强 Ridge 正则化 + 数值稳定方法
  - 评估: 5-fold CV, 不按 group (因为数据已经"散")
"""
import pickle
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import KFold

# 注入 dummy 函数以防 pickle 反序列化报错
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

# ---------- 1. 解析 ----------
print("Parsing all mutations (no WT validation)...")
df["parsed_muts"] = df["aaMutations"].apply(parse_mut_v2)
df["n_mut"] = df["parsed_muts"].apply(len)
# 过滤: 排除 n_mut=0 (WT 单独处理) 和 无法解析的
df_valid = df[df["n_mut"] > 0].copy()
df_valid = df_valid[df_valid["parsed_muts"].apply(len) > 0]
print(f"  Valid mutated rows: {len(df_valid)}")

# ---------- 2. 构建 (type, pos, new_aa) 索引 ----------
print("Building vocab...")
mut_vocab = defaultdict(set)
for _, r in df_valid.iterrows():
    t = r["GFP type"]
    for pos, wt_aa, new_aa in r["parsed_muts"]:
        mut_vocab[t].add((pos, new_aa))

for t in mut_vocab:
    print(f"  {t}: {len(mut_vocab[t])} unique (pos, new_aa)")

# ---------- 3. 拟合模型 (per GFP type) ----------
MAX_NMUT = 8
df_train = df_valid[df_valid["n_mut"] <= MAX_NMUT].copy()
print(f"\nTraining on rows with 1 <= n_mut <= {MAX_NMUT}: {len(df_train)}")

models = {}
for t in ["avGFP", "amacGFP", "cgreGFP", "ppluGFP"]:
    print(f"\n=== {t} ===")
    sub = df_train[df_train["GFP type"] == t].reset_index(drop=True)
    n = len(sub)
    if n < 50:
        print(f"  Skip (n={n} too small)")
        continue

    # 只保留出现 >=3 次的 (pos, new_aa) — 减少噪声
    from collections import Counter
    mut_counter = Counter()
    for pm in sub["parsed_muts"]:
        for k in pm:
            mut_counter[k] += 1
    vocab = sorted([k for k, c in mut_counter.items() if c >= 3])
    vocab_idx = {v: i for i, v in enumerate(vocab)}
    v = len(vocab)
    print(f"  n={n}, vocab(filtered freq>=3)={v}")

    X = np.zeros((n, v), dtype=np.float32)
    y = sub["Brightness"].values
    for i, pm in enumerate(sub["parsed_muts"]):
        for k in pm:
            if k in vocab_idx:
                X[i, vocab_idx[k]] = 1.0

    # 5-fold CV
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    cv_r2, cv_rmse = [], []
    for fold, (tr, va) in enumerate(kf.split(X)):
        m = RidgeCV(alphas=[0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0])
        m.fit(X[tr], y[tr])
        pred = m.predict(X[va])
        cv_r2.append(r2_score(y[va], pred))
        cv_rmse.append(np.sqrt(mean_squared_error(y[va], pred)))
    print(f"  CV R2 = {np.mean(cv_r2):+.3f} (+/- {np.std(cv_r2):.3f})")
    print(f"  CV RMSE = {np.mean(cv_rmse):.3f}")

    # 全数据模型
    m_full = RidgeCV(alphas=[0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0])
    m_full.fit(X, y)
    intercept = float(m_full.intercept_)
    effects = {vocab[i]: float(m_full.coef_[i]) for i in range(v)}

    sorted_eff = sorted(effects.items(), key=lambda x: -x[1])
    print(f"  Top 20 positive effects:")
    for key, e in sorted_eff[:20]:
        if len(key) == 2:
            pos, new_aa = key
        else:
            pos, wt_aa, new_aa = key
        wt_aa_ref = wt_seqs.get(t, '??????')[pos-1] if pos-1 < len(wt_seqs.get(t, '')) else '?'
        print(f"    pos{pos:>3} {wt_aa_ref}->{new_aa}: {e:+.3f}")
    print(f"  Top 10 negative effects:")
    for key, e in sorted_eff[-10:]:
        if len(key) == 2:
            pos, new_aa = key
        else:
            pos, wt_aa, new_aa = key
        wt_aa_ref = wt_seqs.get(t, '??????')[pos-1] if pos-1 < len(wt_seqs.get(t, '')) else '?'
        print(f"    pos{pos:>3} {wt_aa_ref}->{new_aa}: {e:+.3f}")

    models[t] = {
        "intercept": intercept,
        "effects": effects,
        "alpha": float(m_full.alpha_),
        "vocab": vocab,
        "n_train": n,
        "cv_r2_mean": float(np.mean(cv_r2)),
        "cv_rmse_mean": float(np.mean(cv_rmse)),
    }

# ---------- 4. 保存 ----------
with open(PHASE1 / "additive_models_v2.pkl", "wb") as f:
    pickle.dump(models, f)
print(f"\n[OK] Models saved to {PHASE1 / 'additive_models_v2.pkl'}")

# JSON 报告
import json
report = {
    "wt_brightness": wt_brightness,
    "models": {
        t: {
            "intercept": m["intercept"],
            "alpha": m["alpha"],
            "n_train": m["n_train"],
            "cv_r2_mean": m["cv_r2_mean"],
            "cv_rmse_mean": m["cv_rmse_mean"],
            "n_vocab": len(m["vocab"]),
            "top_positive": [
                {"pos": int(k[0]), "new_aa": k[1], "wt_aa": wt_seqs[t][k[0]-1] if k[0]-1 < len(wt_seqs[t]) else "?", "effect": round(v, 4)}
                for k, v in sorted(m["effects"].items(), key=lambda x: -x[1])[:30]
            ],
            "top_negative": [
                {"pos": int(k[0]), "new_aa": k[1], "wt_aa": wt_seqs[t][k[0]-1] if k[0]-1 < len(wt_seqs[t]) else "?", "effect": round(v, 4)}
                for k, v in sorted(m["effects"].items(), key=lambda x: x[1])[:30]
            ],
        } for t, m in models.items()
    }
}
with open(PHASE1 / "model_report_v2.json", "w") as f:
    json.dump(report, f, indent=2)
print(f"[OK] Report saved to {PHASE1 / 'model_report_v2.json'}")
print("\n=== Phase 1.3 v2 DONE ===")
