"""
Round 5 P0-2: EVOLVEpro 风格 RF 回归器
==================================
依据: Science 2024 (Jiang et al.) — PLM (ESM-2) + top-layer RF (Random Forest)
- 用比赛官方 141K 数据训练
- 给 71 条候选独立打分
- 解决 Round 2 OOD 失败问题 (XGBoost val R²=0.92 但对论文突变低估)

策略:
  1. 用 Round 2 已有 esm650m_embeddings.npy (141K × 1280)
  2. 训练 RF 回归器 (EVOLVEpro 论文最优)
  3. 对 71 个候选 (Round 4 v5 + Round 5 LMPNN) 计算 ESM-2 嵌入
  4. RF 预测亮度
"""
import numpy as np
import pandas as pd
import json
import time
from pathlib import Path
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score

ROOT = Path(r"D:\生信\2026Protein Design")
R2 = ROOT / "work" / "round2"
R5 = ROOT / "work" / "round5"

EMB_PATH = R2 / "esm650m_embeddings.npy"
IDS_PATH = R2 / "esm650m_ids.csv"

print("=" * 70)
print("EVOLVEpro 风格 fine-tune (ESM-2 嵌入 + RF top layer)")
print("=" * 70)

# 1. 加载训练数据
print("\n1. 加载训练数据...")
emb = np.load(EMB_PATH, mmap_mode="r")
ids = pd.read_csv(IDS_PATH)
print(f"  embeddings shape={emb.shape}")
print(f"  ids cols={list(ids.columns)}")
print(f"  brightness range: [{ids['brightness'].min():.3f}, {ids['brightness'].max():.3f}]")

types = ids["type"].astype(str).values
unique_types = sorted(set(types))
type_to_idx = {t: i for i, t in enumerate(unique_types)}
type_oh = np.eye(len(unique_types), dtype=np.float32)[np.array([type_to_idx[t] for t in types])]

# 用纯 ESM 嵌入 + type one-hot (与 Round 2 一致)
X = np.concatenate([emb[:].astype(np.float32), type_oh], axis=1)
y = ids["brightness"].astype(np.float32).values
print(f"  X shape={X.shape}, y shape={y.shape}")

# 2. 训练 RF (EVOLVEpro 论文最优配置)
# 论文超参: n_estimators=300, max_depth=None, min_samples_split=2
print("\n2. 训练 Random Forest 回归器...")
t0 = time.time()
rf = RandomForestRegressor(
    n_estimators=200,  # 减少以加速
    max_depth=20,
    min_samples_split=5,
    n_jobs=-1,
    random_state=42,
    verbose=0,
)

# 9:1 train/val split (与 Round 2 一致)
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

print("  训练中... (预计 5-10 min)")
rf.fit(X_train, y_train)

# 评估
y_pred_train = rf.predict(X_train)
y_pred_val = rf.predict(X_val)
train_r2 = r2_score(y_train, y_pred_train)
val_r2 = r2_score(y_val, y_pred_val)
elapsed = time.time() - t0
print(f"  ✓ 训练完成 ({elapsed:.1f}s)")
print(f"  train R²: {train_r2:.4f}")
print(f"  val R²:   {val_r2:.4f}")

# 保存模型
import pickle
with open(R5 / "evolvepro_rf.pkl", "wb") as f:
    pickle.dump({"rf": rf, "type_to_idx": type_to_idx, "val_r2": val_r2, "train_r2": train_r2}, f)
print(f"  模型保存到 evolvepro_rf.pkl")

# 3. 给候选评分
print("\n3. 候选 ESM-2 嵌入 + RF 预测...")
import torch
from transformers import AutoTokenizer, AutoModel

# 加载所有 Round 4 + Round 5 候选
with open(R5 / "final_6_round5.json", encoding="utf-8") as f:
    final_6 = json.load(f)

# 也加载所有候选 (71 条) 做完整评分
with open(ROOT / "work" / "round4" / "esmfold_round4_v3.json", encoding="utf-8") as f:
    r4_hc = json.load(f)
with open(ROOT / "work" / "round4" / "esmfold_mpnn.json", encoding="utf-8") as f:
    r4_msf = json.load(f)
with open(ROOT / "work" / "round4" / "esmfold_mpnn_av.json", encoding="utf-8") as f:
    r4_mav = json.load(f)
with open(R5 / "esmfold_lmpnn_v2.json", encoding="utf-8") as f:
    r5_lm = json.load(f)

all_cands = r4_hc + r4_msf + r4_mav + r5_lm
print(f"  总候选: {len(all_cands)}")

# 提取 ESM-2 650M 嵌入
print("  加载 ESM-2 650M...")
tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t33_650M_UR50D", local_files_only=True)
model = AutoModel.from_pretrained("facebook/esm2_t33_650M_UR50D", local_files_only=True).cuda().eval()

def embed_seq(seq):
    with torch.no_grad():
        inputs = tokenizer(seq, return_tensors="pt", add_special_tokens=False)
        inputs = {k: v.cuda() for k, v in inputs.items()}
        out = model(**inputs)
        # 全序列平均池化 (与训练时一致)
        return out.last_hidden_state[0].mean(dim=0).cpu().numpy()

# 给候选添加 type one-hot (匹配训练特征)
# 候选骨架映射到训练时的类型
scaffold_to_type = {
    "sfGFP": "avGFP",  # sfGFP 是 avGFP 的衍生
    "avGFP": "avGFP",
    "amacGFP": "amacGFP",
    "cgreGFP": "cgreGFP",
    "ppluGFP": "ppluGFP",
    "mBaoJin": "avGFP",  # 假设映射到 avGFP (训练集无 StayGold)
    "sfGFP_MPNN": "avGFP",
    "avGFP_MPNN": "avGFP",
    "avGFP_LMPNN": "avGFP",
}

# 批量推理
print(f"  对 {len(all_cands)} 个候选做嵌入推理...")
t0 = time.time()
for i, c in enumerate(all_cands):
    if i % 10 == 0:
        print(f"    {i}/{len(all_cands)}...")
    seq = c["seq"]
    emb_vec = embed_seq(seq)
    
    type_str = scaffold_to_type.get(c["scaffold"], "avGFP")
    type_idx = type_to_idx.get(type_str, 0)
    type_vec = np.zeros(len(unique_types), dtype=np.float32)
    type_vec[type_idx] = 1.0
    
    feat = np.concatenate([emb_vec, type_vec])[None, :]
    pred = rf.predict(feat)[0]
    c["evolvepro_pred_brightness"] = float(pred)

print(f"  完成 ({time.time()-t0:.1f}s)")

# 排序展示
print("\n4. Top-15 候选 (按 EVOLVEpro 预测亮度):")
print("=" * 100)
sorted_cands = sorted(all_cands, key=lambda x: -x.get("evolvepro_pred_brightness", 0))
print(f"{'#':<3} {'name':<30} {'scaffold':<14} {'pLDDT':>5} {'EVOLVEpro_brt':>12}")
print("-" * 100)
for i, c in enumerate(sorted_cands[:15], 1):
    print(f"{i:<3} {c['name'][:30]:<30} {c['scaffold'][:14]:<14} "
          f"{c['plddt_mean']:>5.1f} {c['evolvepro_pred_brightness']:>12.3f}")

# WT 参考
print(f"\n参考 WT 训练集 log10 brightness:")
for t in unique_types:
    t_mask = types == t
    print(f"  {t}: mean={y[t_mask].mean():.3f}, max={y[t_mask].max():.3f}")
print(f"  全体: mean={y.mean():.3f}, max={y.max():.3f}")

with open(R5 / "evolvepro_scored.json", "w", encoding="utf-8") as f:
    json.dump(sorted_cands, f, indent=2, ensure_ascii=False)
print(f"\n✓ 保存 evolvepro_scored.json")
