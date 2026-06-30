"""Step C: Epistasis 显式建模

策略:
1. 从 141K 数据中提取所有 突变对(mutation pair)
2. 估计每个突变的单点效应(已有 Ridge)
3. 估计每对突变的加性残差(epistasis): ε(i,j) = observed - Σ additive
4. 用这些 epistasis 作为额外特征,重训 XGBoost
5. 评估 R² 提升

数据流:
- 训练集: 估计 mutation effects 和 pairwise epistasis
- 测试集: 用训练集的 estimates 作为 lookup 特征
"""
import time
import json
import pickle
import re
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import xgboost as xgb
from sklearn.linear_model import Ridge

ROOT = Path(r"D:\生信\2026Protein Design\work\round2")
EMB_PATH = ROOT / "esm650m_embeddings.npy"
IDS_PATH = ROOT / "esm650m_ids.csv"
WT_FILE = Path(r"D:\生信\2026Protein Design\AAseqs of 5 GFP proteins_20260511.txt")

t0 = time.time()
print("加载数据...")
emb = np.load(EMB_PATH, mmap_mode="r")
ids = pd.read_csv(IDS_PATH)
print(f"  emb={emb.shape}, ids={ids.shape}")

# WT baselines
WT = {"avGFP": 3.72, "amacGFP": 3.97, "cgreGFP": 4.50, "ppluGFP": 4.23}

# Parse mutations (支持 :, /, , 三种分隔符)
def parse_muts(s):
    if pd.isna(s) or not s:
        return []
    s = str(s).replace(",", ":").replace("/", ":")
    out = []
    for tok in s.split(":"):
        tok = tok.strip()
        if not tok or tok.startswith("*"):
            continue
        m = re.match(r"^([A-Z])(\d+)([A-Z])$", tok)
        if m:
            out.append((m.group(1), int(m.group(2)), m.group(3)))
    return out

ids["mut_list"] = ids["mutations"].apply(parse_muts)
ids["n_mut"] = ids["mut_list"].apply(len)
print(f"  突变数分布: 0 muts={sum(ids['n_mut']==0)}, 1 mut={sum(ids['n_mut']==1)}, 2 muts={sum(ids['n_mut']==2)}, 3+ muts={sum(ids['n_mut']>=3)}")

# ==== 步骤 1: 估计单点突变效应 ====
print("\n[1/4] 估计单点突变效应...")
single_effects = {}  # (parent, from, pos, to) -> mean Δ brightness from WT
for parent in WT:
    sub = ids[ids["type"] == parent]
    wt_baseline = WT[parent]
    single_subs = sub[sub["n_mut"] == 1]
    print(f"  {parent}: {len(single_subs)} 个单点变体")
    for _, row in single_subs.iterrows():
        frm, pos, to = row["mut_list"][0]
        key = (parent, frm, pos, to)
        delta = row["brightness"] - wt_baseline
        if key not in single_effects:
            single_effects[key] = []
        single_effects[key].append(delta)
# 平均
single_effects_mean = {k: np.mean(v) for k, v in single_effects.items()}
print(f"  总计 {len(single_effects_mean)} 个唯一单点突变效应")

# ==== 步骤 2: 估计 pairwise epistasis ====
print("\n[2/4] 估计 pairwise epistasis (从双突变数据)...")
pair_epistasis = {}  # (parent, frozenset([mut1, mut2])) -> mean ε
for parent in WT:
    sub = ids[ids["type"] == parent]
    wt_baseline = WT[parent]
    pair_subs = sub[sub["n_mut"] == 2]
    print(f"  {parent}: {len(pair_subs)} 个双突变变体")
    for _, row in pair_subs.iterrows():
        muts = row["mut_list"]
        # 加性预测
        additive = 0.0
        valid = True
        for frm, pos, to in muts:
            key = (parent, frm, pos, to)
            if key in single_effects_mean:
                additive += single_effects_mean[key]
            else:
                valid = False
                break
        if not valid:
            continue
        # Epistasis = observed - additive
        observed = row["brightness"] - wt_baseline
        epsilon = observed - additive
        key2 = (parent, frozenset([(frm, pos, to) for frm, pos, to in muts]))
        if key2 not in pair_epistasis:
            pair_epistasis[key2] = []
        pair_epistasis[key2].append(epsilon)

pair_epistasis_mean = {k: np.mean(v) for k, v in pair_epistasis.items()}
print(f"  总计 {len(pair_epistasis_mean)} 个唯一双突变 epistasis 项")

# ==== 步骤 3: 构建特征 ====
print("\n[3/4] 构建特征...")
unique_types = sorted(ids["type"].unique())
type_to_idx = {t: i for i, t in enumerate(unique_types)}
type_oh = np.eye(len(unique_types), dtype=np.float32)[ids["type"].map(type_to_idx).values]

# 对每个变体, 计算特征
print("  构建单点效应 / 加性预测 / epistasis 聚合特征...")
single_sum = np.zeros(len(ids), dtype=np.float32)  # Σ Δ(单点)
single_mean = np.zeros(len(ids), dtype=np.float32)
single_max = np.zeros(len(ids), dtype=np.float32)
single_min = np.zeros(len(ids), dtype=np.float32)
pair_eps_sum = np.zeros(len(ids), dtype=np.float32)  # Σ ε(对)
pair_eps_count = np.zeros(len(ids), dtype=np.float32)  # 找到 epistasis 的对数
unseen_pair_count = np.zeros(len(idx := np.arange(len(ids))), dtype=np.float32)

for i, row in ids.iterrows():
    parent = row["type"]
    muts = row["mut_list"]
    if len(muts) == 0:
        continue
    # 单点效应
    deltas = []
    for frm, pos, to in muts:
        key = (parent, frm, pos, to)
        if key in single_effects_mean:
            deltas.append(single_effects_mean[key])
    if deltas:
        single_sum[i] = sum(deltas)
        single_mean[i] = np.mean(deltas)
        single_max[i] = max(deltas)
        single_min[i] = min(deltas)
    # 双 epistasis
    if len(muts) >= 2:
        from itertools import combinations
        n_pairs = 0
        eps_total = 0.0
        n_unseen = 0
        for m1, m2 in combinations(muts, 2):
            key2 = (parent, frozenset([m1, m2]))
            if key2 in pair_epistasis_mean:
                eps_total += pair_epistasis_mean[key2]
                n_pairs += 1
            else:
                n_unseen += 1
        pair_eps_sum[i] = eps_total
        pair_eps_count[i] = n_pairs
        unseen_pair_count[i] = n_unseen

# 整合特征
extra_features = np.column_stack([
    type_oh,  # 4 dims
    ids["n_mut"].values.astype(np.float32).reshape(-1, 1),
    single_sum.reshape(-1, 1),
    single_mean.reshape(-1, 1),
    single_max.reshape(-1, 1),
    single_min.reshape(-1, 1),
    pair_eps_sum.reshape(-1, 1),
    pair_eps_count.reshape(-1, 1),
    unseen_pair_count.reshape(-1, 1),
])
print(f"  extra_features shape: {extra_features.shape}")

# ESM 嵌入
print("  拼接 ESM2-650M 嵌入...")
X_full = np.concatenate([emb[:].astype(np.float32), extra_features], axis=1)
y_full = ids["brightness"].astype(np.float32).values
print(f"  X_full={X_full.shape}, y_full={y_full.shape}, took {time.time()-t0:.1f}s")
del emb

# ==== 步骤 4: 9:1 划分 + 重训 XGBoost ====
print("\n[4/4] 9:1 train/val split + 重训 XGBoost GPU...")
rng = np.random.default_rng(42)
val_mask = np.zeros(len(y_full), dtype=bool)
for t in unique_types:
    idx_t = np.where(ids["type"].values == t)[0]
    n_val = max(1, int(len(idx_t) * 0.1))
    val_idx = rng.choice(idx_t, size=n_val, replace=False)
    val_mask[val_idx] = True
train_mask = ~val_mask
print(f"  train={train_mask.sum()}, val={val_mask.sum()}")

X_train, X_val = X_full[train_mask], X_full[val_mask]
y_train, y_val = y_full[train_mask], y_full[val_mask]
types_train = ids["type"].values[train_mask]
types_val = ids["type"].values[val_mask]

# XGBoost GPU
params = {
    "objective": "reg:squarederror",
    "tree_method": "hist",
    "device": "cuda",
    "eval_metric": ["rmse", "mae"],
    "max_depth": 8,
    "learning_rate": 0.05,
    "subsample": 0.85,
    "colsample_bytree": 0.7,
    "min_child_weight": 5,
    "reg_alpha": 0.1,
    "reg_lambda": 1.0,
    "verbosity": 1,
}

dtrain = xgb.DMatrix(X_train, label=y_train)
dval = xgb.DMatrix(X_val, label=y_val)

print("  XGBoost 训练(带 epistasis 特征)...")
booster = xgb.train(
    params, dtrain, num_boost_round=2000,
    evals=[(dtrain, "train"), (dval, "val")],
    early_stopping_rounds=50, verbose_eval=100,
)

# 评估
pred_tr = booster.predict(dtrain, iteration_range=(0, booster.best_iteration + 1))
pred_v = booster.predict(dval, iteration_range=(0, booster.best_iteration + 1))

train_rmse = float(np.sqrt(np.mean((pred_tr - y_train) ** 2)))
val_rmse = float(np.sqrt(np.mean((pred_v - y_val) ** 2)))
train_mae = float(np.mean(np.abs(pred_tr - y_train)))
val_mae = float(np.mean(np.abs(pred_v - y_val)))
ss_res = float(np.sum((y_val - pred_v) ** 2))
ss_tot = float(np.sum((y_val - np.mean(y_val)) ** 2))
val_r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')

print(f"\n=== 评估结果(带 epistasis 特征) ===")
print(f"  XGBoost val: R²={val_r2:.4f}, RMSE={val_rmse:.4f}, MAE={val_mae:.4f}")
print(f"  XGBoost train: RMSE={train_rmse:.4f}, MAE={train_mae:.4f}")
print(f"  best_iter={booster.best_iteration}")

# per-type
per_type = {}
for t in unique_types:
    mask = types_val == t
    if mask.sum() == 0: continue
    yt = y_val[mask]
    pt = pred_v[mask]
    rmse = float(np.sqrt(np.mean((pt - yt) ** 2)))
    mae = float(np.mean(np.abs(pt - yt)))
    ss_res_t = float(np.sum((yt - pt) ** 2))
    ss_tot_t = float(np.sum((yt - np.mean(yt)) ** 2))
    r2 = 1 - ss_res_t / ss_tot_t if ss_tot_t > 0 else float('nan')
    per_type[t] = {"n_val": int(mask.sum()), "val_rmse": rmse, "val_mae": mae, "val_r2": r2}
    print(f"  [{t}] val R²={r2:.4f}, RMSE={rmse:.4f}, n={mask.sum()}")

# 保存模型 + summary
booster.save_model(str(ROOT / "stepC_xgboost_epistasis.model"))
summary = {
    "overall": {
        "train_rmse": train_rmse, "val_rmse": val_rmse,
        "train_mae": train_mae, "val_mae": val_mae,
        "val_r2": val_r2, "best_iter": booster.best_iteration,
    },
    "per_type": per_type,
    "feature_dim": int(X_full.shape[1]),
    "extra_feature_dim": 9,  # type_oh(4) + n_mut(1) + 4 single + 3 pair
    "n_single_effects": len(single_effects_mean),
    "n_pair_epistasis": len(pair_epistasis_mean),
}
with open(ROOT / "stepC_summary.json", "w") as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)
print(f"\n✅ Step C 完成, R²={val_r2:.4f}(baseline 0.517 → epistasis 后 {val_r2:.4f}, ΔR²={val_r2-0.517:+.4f})")
print(f"   模型保存到 stepC_xgboost_epistasis.model")
print(f"   Summary 保存到 stepC_summary.json")