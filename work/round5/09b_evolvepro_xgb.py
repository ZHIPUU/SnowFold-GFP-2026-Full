"""
Round 5 P0-2 优化版: EVOLVEpro 风格 fine-tune (XGBoost GPU)
==================================
为什么改用 XGBoost GPU:
  - RandomForest 仅 CPU, 140K × 1280 维特征训练 > 10 min
  - XGBoost GPU 同样作为 top-layer 回归器 (EVOLVEpro 论文允许),
    速度快 10x, Round 2 已验证可行
  - 目标: 给 98 候选独立打分, 解决 OOD 问题

策略 (与 Round 2 baseline R²=0.517 区分):
  - 使用 GPU 加速
  - 优化超参 (max_depth=6, lr=0.1, 减小过拟合)
  - 但避免 Round 2 epistasis 特征 (那才是 R²=0.916 但 OOD 失败的元凶)
"""
import numpy as np
import pandas as pd
import json
import time
import pickle
from pathlib import Path
from sklearn.metrics import r2_score
import xgboost as xgb

ROOT = Path(r"D:\生信\2026Protein Design")
R2 = ROOT / "work" / "round2"
R5 = ROOT / "work" / "round5"

print("=" * 70)
print("EVOLVEpro 风格 (ESM-2 嵌入 + XGBoost GPU top layer)")
print("=" * 70)

# 1. 数据
print("\n1. 加载数据...")
emb = np.load(R2 / "esm650m_embeddings.npy", mmap_mode="r")
ids = pd.read_csv(R2 / "esm650m_ids.csv")
print(f"  embeddings: {emb.shape}")

types = ids["type"].astype(str).values
unique_types = sorted(set(types))
type_to_idx = {t: i for i, t in enumerate(unique_types)}
type_oh = np.eye(len(unique_types), dtype=np.float32)[np.array([type_to_idx[t] for t in types])]

X = np.concatenate([emb[:].astype(np.float32), type_oh], axis=1)
y = ids["brightness"].astype(np.float32).values
print(f"  X: {X.shape}, y range: [{y.min():.3f}, {y.max():.3f}]")

# 2. Train/val split
print("\n2. 9:1 train/val split (stratified)...")
rng = np.random.default_rng(42)
val_mask = np.zeros(len(y), dtype=bool)
for t in unique_types:
    type_idx = np.where(types == t)[0]
    n_val = len(type_idx) // 10
    val_idx = rng.choice(type_idx, size=n_val, replace=False)
    val_mask[val_idx] = True

X_train, y_train = X[~val_mask], y[~val_mask]
X_val, y_val = X[val_mask], y[val_mask]
print(f"  train: {len(X_train)}, val: {len(X_val)}")

# 3. XGBoost GPU 训练
print("\n3. 训练 XGBoost GPU...")
t0 = time.time()
params = {
    "objective": "reg:squarederror",
    "tree_method": "hist",
    "device": "cuda",
    "max_depth": 6,            # 浅一些, 减少过拟合
    "learning_rate": 0.1,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "min_child_weight": 5,
    "reg_alpha": 0.1,
    "reg_lambda": 1.0,
    "n_estimators": 500,
}

dtrain = xgb.DMatrix(X_train, label=y_train)
dval = xgb.DMatrix(X_val, label=y_val)

booster = xgb.train(
    params, dtrain,
    num_boost_round=500,
    evals=[(dtrain, "train"), (dval, "val")],
    early_stopping_rounds=30,
    verbose_eval=50,
)
elapsed = time.time() - t0
print(f"  训练完成: {elapsed:.1f}s, best iter: {booster.best_iteration}")

y_pred_val = booster.predict(dval)
y_pred_train = booster.predict(dtrain)
train_r2 = r2_score(y_train, y_pred_train)
val_r2 = r2_score(y_val, y_pred_val)
print(f"  train R²: {train_r2:.4f}")
print(f"  val R²:   {val_r2:.4f}")

booster.save_model(str(R5 / "evolvepro_xgb.model"))
print(f"  模型保存到 evolvepro_xgb.model")

# 4. 加载所有 98 候选 + ESM-2 嵌入
print("\n4. 候选 ESM-2 嵌入...")
import torch
from transformers import AutoTokenizer, AutoModel

# 加载所有候选
with open(R5 / "final_6_round5_v2.json", encoding="utf-8") as f:
    final_6 = json.load(f)

with open(ROOT / "work" / "round4" / "esmfold_round4_v3.json", encoding="utf-8") as f:
    r4_hc = json.load(f)
with open(ROOT / "work" / "round4" / "esmfold_mpnn.json", encoding="utf-8") as f:
    r4_msf = json.load(f)
with open(ROOT / "work" / "round4" / "esmfold_mpnn_av.json", encoding="utf-8") as f:
    r4_mav = json.load(f)
with open(R5 / "esmfold_lmpnn_v2.json", encoding="utf-8") as f:
    r5_lm = json.load(f)
with open(R5 / "esmfold_lmpnn_expanded.json", encoding="utf-8") as f:
    r5_le = json.load(f)

all_cands = r4_hc + r4_msf + r4_mav + r5_lm + r5_le
# 去重
seen = set()
dedup = []
for c in all_cands:
    if c["seq"] not in seen:
        seen.add(c["seq"])
        dedup.append(c)
print(f"  总候选 (去重): {len(dedup)}")

# ESM-2 650M 嵌入
print("  加载 ESM-2 650M...")
tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t33_650M_UR50D", local_files_only=True)
model = AutoModel.from_pretrained("facebook/esm2_t33_650M_UR50D", local_files_only=True).cuda().eval()
print("  模型就绪 (GPU)")

scaffold_to_type = {
    "sfGFP": "avGFP", "avGFP": "avGFP", "amacGFP": "amacGFP",
    "cgreGFP": "cgreGFP", "ppluGFP": "ppluGFP",
    "mBaoJin": "avGFP",
    "sfGFP_MPNN": "avGFP", "avGFP_MPNN": "avGFP", "avGFP_LMPNN": "avGFP",
}

print(f"  对 {len(dedup)} 个候选做嵌入推理 (GPU)...")
t0 = time.time()
embs_pred = []
for i, c in enumerate(dedup):
    if i % 20 == 0:
        print(f"    {i}/{len(dedup)}... ({time.time()-t0:.1f}s)")
    seq = c["seq"]
    with torch.no_grad():
        inputs = tokenizer(seq, return_tensors="pt", add_special_tokens=False)
        inputs = {k: v.cuda() for k, v in inputs.items()}
        out = model(**inputs)
        emb_vec = out.last_hidden_state[0].mean(dim=0).cpu().numpy()
    
    t_str = scaffold_to_type.get(c["scaffold"], "avGFP")
    t_idx = type_to_idx.get(t_str, 0)
    t_vec = np.zeros(len(unique_types), dtype=np.float32)
    t_vec[t_idx] = 1.0
    
    feat = np.concatenate([emb_vec, t_vec])[None, :]
    pred = booster.predict(xgb.DMatrix(feat))[0]
    c["evolvepro_pred"] = float(pred)

print(f"  完成: {time.time()-t0:.1f}s")

# 5. 排序展示
print("\n5. Top-15 (按 EVOLVEpro 预测亮度):")
print("=" * 100)
sorted_cands = sorted(dedup, key=lambda x: -x.get("evolvepro_pred", 0))
print(f"{'#':<3} {'name':<32} {'scaffold':<14} {'pLDDT':>5} {'EVOLVEpro':>9}")
print("-" * 80)
for i, c in enumerate(sorted_cands[:15], 1):
    print(f"{i:<3} {c['name'][:32]:<32} {c['scaffold'][:14]:<14} "
          f"{c['plddt_mean']:>5.1f} {c['evolvepro_pred']:>9.3f}")

# WT 参考
print(f"\nWT 训练集参考 (log10 brightness):")
for t in unique_types:
    t_mask = types == t
    print(f"  {t}: mean={y[t_mask].mean():.3f}, max={y[t_mask].max():.3f}")

with open(R5 / "evolvepro_scored.json", "w", encoding="utf-8") as f:
    json.dump(sorted_cands, f, indent=2, ensure_ascii=False)
print(f"\n✓ 保存 evolvepro_scored.json")
print(f"  val_R² = {val_r2:.4f}")
