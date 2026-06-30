"""Step 2: 用 ESM2-650M 嵌入训练亮度预测模型
- 主模型: XGBoost GPU (1280-dim 嵌入 + type one-hot)
- Baseline: 加性 Ridge (4 个母体单独,纯突变加性)
- 输出 R² / MAE / RMSE (train/val)
"""
import time
import numpy as np
import pandas as pd
from pathlib import Path
import pickle
import json

ROOT = Path(r"D:\生信\2026Protein Design\work\round2")
EMB_PATH = ROOT / "esm650m_embeddings.npy"
IDS_PATH = ROOT / "esm650m_ids.csv"
LOG_PATH = ROOT / "step2_log.txt"

def log(msg):
    s = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(s, flush=True)
    with open(LOG_PATH, "a", encoding="utf-8") as f:
        f.write(s + "\n")

log("=== Step 2: 训练 XGBoost GPU + 加性 Ridge baseline ===")

# 1) 加载数据
log("加载嵌入 + 元数据...")
emb = np.load(EMB_PATH, mmap_mode="r")
ids = pd.read_csv(IDS_PATH)
log(f"  embeddings shape={emb.shape}, dtype={emb.dtype}")
log(f"  ids shape={ids.shape}, cols={list(ids.columns)}")

# 2) 准备特征 X 和目标 y
log("构建特征矩阵 X (嵌入 + type one-hot)...")
types = ids["type"].astype(str).values
unique_types = sorted(set(types))
type_to_idx = {t: i for i, t in enumerate(unique_types)}
type_oh = np.eye(len(unique_types), dtype=np.float32)[np.array([type_to_idx[t] for t in types])]
log(f"  types={unique_types}, type_oh shape={type_oh.shape}")

# 用 ESM 嵌入 (1280-d) + type one-hot (4-d)
X = np.concatenate([emb[:].astype(np.float32), type_oh], axis=1)
y = ids["brightness"].astype(np.float32).values  # log10 单位
log(f"  X shape={X.shape}, y shape={y.shape}, y range=[{y.min():.3f},{y.max():.3f}]")

# 3) 9:1 train/val split (stratified by type)
log("9:1 train/val split (stratified by type)...")
rng = np.random.default_rng(42)
val_mask = np.zeros(len(y), dtype=bool)
for t in unique_types:
    idx = np.where(types == t)[0]
    n_val = max(1, int(len(idx) * 0.1))
    val_idx = rng.choice(idx, size=n_val, replace=False)
    val_mask[val_idx] = True
train_mask = ~val_mask
log(f"  train={train_mask.sum()}, val={val_mask.sum()}")

X_train, X_val = X[train_mask], X[val_mask]
y_train, y_val = y[train_mask], y[val_mask]
types_train, types_val = types[train_mask], types[val_mask]
log(f"  X_train={X_train.shape}, X_val={X_val.shape}")

# 释放 memmap
del emb, X

# 4) 加性 Ridge baseline (per-parent, 突变加性)
log("\n--- Baseline: 加性 Ridge (per-parent) ---")
import re
from sklearn.linear_model import Ridge

# 解析 mutations 字符串
def parse_muts(s):
    """支持 : / , 三种分隔符"""
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

mutations = ids["mutations"].apply(parse_muts).values
log(f"  平均突变数/变体: {np.mean([len(m) for m in mutations]):.2f}")

# 对每个 type 训练 Ridge
baseline_preds_val = np.full_like(y_val, np.nan)
baseline_preds_train = np.full_like(y_train, np.nan)
type_perf = {}
for t in unique_types:
    train_t = (types_train == t)
    val_t = (types_val == t)
    if train_t.sum() < 50:
        log(f"  {t}: 跳过 (train n={train_t.sum()})")
        continue

    # 收集所有突变位置
    mut_idx_train = np.where(train_t)[0]
    mut_idx_val = np.where(val_t)[0]
    all_muts = set()
    for mi in mut_idx_train:
        for aa_from, pos, aa_to in mutations[mi]:
            all_muts.add((pos, aa_from, aa_to))
    all_muts = sorted(all_muts)
    log(f"  {t}: train={train_t.sum()}, val={val_t.sum()}, 唯一突变数={len(all_muts)}")

    # 特征矩阵: 0/1 是否含此突变
    mut_to_col = {m: i for i, m in enumerate(all_muts)}
    n_feat = len(all_muts)
    Xt = np.zeros((train_t.sum(), n_feat), dtype=np.float32)
    Xv = np.zeros((val_t.sum(), n_feat), dtype=np.float32)
    yt = y_train[train_t]
    yv = y_val[val_t]

    for ri, mi in enumerate(mut_idx_train):
        for aa_from, pos, aa_to in mutations[mi]:
            if (pos, aa_from, aa_to) in mut_to_col:
                Xt[ri, mut_to_col[(pos, aa_from, aa_to)]] = 1.0
    for ri, mi in enumerate(mut_idx_val):
        for aa_from, pos, aa_to in mutations[mi]:
            if (pos, aa_from, aa_to) in mut_to_col:
                Xv[ri, mut_to_col[(pos, aa_from, aa_to)]] = 1.0

    # 训练 Ridge
    ridge = Ridge(alpha=1.0)
    ridge.fit(Xt, yt)
    pred_tr = ridge.predict(Xt)
    pred_v = ridge.predict(Xv)

    # 反 log10: 看相对亮度
    rel_tr = np.mean(10 ** pred_tr) / np.mean(10 ** yt)
    rel_v = np.mean(10 ** pred_v) / np.mean(10 ** yv)

    train_rmse = float(np.sqrt(np.mean((pred_tr - yt) ** 2)))
    val_rmse = float(np.sqrt(np.mean((pred_v - yv) ** 2)))
    train_mae = float(np.mean(np.abs(pred_tr - yt)))
    val_mae = float(np.mean(np.abs(pred_v - yv)))

    ss_res = float(np.sum((yv - pred_v) ** 2))
    ss_tot = float(np.sum((yv - np.mean(yv)) ** 2))
    val_r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')

    type_perf[t] = {
        "n_train": int(train_t.sum()),
        "n_val": int(val_t.sum()),
        "n_muts": n_feat,
        "train_rmse": train_rmse,
        "val_rmse": val_rmse,
        "train_mae": train_mae,
        "val_mae": val_mae,
        "val_r2": val_r2,
        "val_rel_brightness": float(rel_v),
    }
    log(f"    Ridge val R²={val_r2:.3f}, RMSE={val_rmse:.3f}, MAE={val_mae:.3f}, rel_brightness={rel_v:.2f}×")

    # 填回大数组
    train_t_idx = np.where(train_t)[0]
    val_t_idx = np.where(val_t)[0]
    baseline_preds_train[train_t_idx] = pred_tr
    baseline_preds_val[val_t_idx] = pred_v

# 总体 baseline
val_rmse_all = float(np.sqrt(np.mean((baseline_preds_val - y_val) ** 2)))
val_mae_all = float(np.mean(np.abs(baseline_preds_val - y_val)))
ss_res = float(np.sum((y_val - baseline_preds_val) ** 2))
ss_tot = float(np.sum((y_val - np.mean(y_val)) ** 2))
val_r2_all = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')
log(f"  Baseline 总体 val: R²={val_r2_all:.3f}, RMSE={val_rmse_all:.3f}, MAE={val_mae_all:.3f}")

baseline_summary = {
    "overall": {"val_r2": val_r2_all, "val_rmse": val_rmse_all, "val_mae": val_mae_all},
    "per_type": type_perf,
}
with open(ROOT / "step2_baseline_ridge.json", "w", encoding="utf-8") as f:
    json.dump(baseline_summary, f, indent=2, ensure_ascii=False)
log(f"  保存 step2_baseline_ridge.json")

# 5) XGBoost GPU (主模型)
log("\n--- 主模型: XGBoost GPU (嵌入 + type one-hot) ---")
import xgboost as xgb

t0 = time.time()
# GPU params
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
log(f"  DMatrix 构建完成 (X_train={X_train.shape}, X_val={X_val.shape})")

# 训练,带早停
booster = xgb.train(
    params,
    dtrain,
    num_boost_round=2000,
    evals=[(dtrain, "train"), (dval, "val")],
    early_stopping_rounds=50,
    verbose_eval=50,
)
t1 = time.time()
log(f"  XGBoost 训练完成,耗时 {t1-t0:.1f}s, best_iter={booster.best_iteration}")

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

# per-type 评估
xgb_per_type = {}
for t in unique_types:
    mask = types_val == t
    if mask.sum() == 0:
        continue
    yt = y_val[mask]
    pt = pred_v[mask]
    rmse = float(np.sqrt(np.mean((pt - yt) ** 2)))
    mae = float(np.mean(np.abs(pt - yt)))
    ss_res_t = float(np.sum((yt - pt) ** 2))
    ss_tot_t = float(np.sum((yt - np.mean(yt)) ** 2))
    r2 = 1 - ss_res_t / ss_tot_t if ss_tot_t > 0 else float('nan')
    rel = float(np.mean(10 ** pt) / np.mean(10 ** yt))
    xgb_per_type[t] = {
        "n_val": int(mask.sum()),
        "val_rmse": rmse,
        "val_mae": mae,
        "val_r2": r2,
        "val_rel_brightness": rel,
    }

xgb_summary = {
    "overall": {
        "train_rmse": train_rmse,
        "val_rmse": val_rmse,
        "train_mae": train_mae,
        "val_mae": val_mae,
        "val_r2": val_r2,
        "best_iter": booster.best_iteration,
        "train_seconds": t1 - t0,
    },
    "per_type": xgb_per_type,
    "params": {k: v for k, v in params.items() if not isinstance(v, list)},
}

log(f"  XGBoost val: R²={val_r2:.3f}, RMSE={val_rmse:.3f}, MAE={val_mae:.3f}")
log(f"  Baseline val: R²={val_r2_all:.3f}, RMSE={val_rmse_all:.3f}, MAE={val_mae_all:.3f}")
log(f"  XGBoost 提升: ΔR²={val_r2 - val_r2_all:+.3f}, RMSE比={val_rmse/val_rmse_all:.2f}")

for t in unique_types:
    if t in xgb_per_type:
        b = baseline_summary["per_type"].get(t, {})
        x = xgb_per_type[t]
        log(f"  [{t}] XGB R²={x['val_r2']:.3f} RMSE={x['val_rmse']:.3f} | Ridge R²={b.get('val_r2',float('nan')):.3f} RMSE={b.get('val_rmse',float('nan')):.3f}")

with open(ROOT / "step2_xgboost_gpu.json", "w", encoding="utf-8") as f:
    json.dump(xgb_summary, f, indent=2, ensure_ascii=False)
log(f"  保存 step2_xgboost_gpu.json")

# 保存 booster
booster.save_model(str(ROOT / "step2_xgboost_gpu.model"))
log(f"  保存 step2_xgboost_gpu.model")

log("\n=== Step 2 DONE ===")
log(f"XGBoost GPU 训练完成: val R²={val_r2:.3f}, RMSE={val_rmse:.3f}, best_iter={booster.best_iteration}")
log(f"Baseline Ridge: val R²={val_r2_all:.3f}, RMSE={val_rmse_all:.3f}")