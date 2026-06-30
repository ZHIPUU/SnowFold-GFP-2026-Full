"""Step 3a: 用新 XGBoost GPU 模型对已有 141K 嵌入打分,找出数据里最强的变体。"""
import time
import json
import numpy as np
import pandas as pd
import xgboost as xgb
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design\work\round2")
EMB = ROOT / "esm650m_embeddings.npy"
IDS = ROOT / "esm650m_ids.csv"
MODEL = ROOT / "step2_xgboost_gpu.model"

t0 = time.time()
print("加载数据...")
emb = np.load(EMB, mmap_mode="r")
ids = pd.read_csv(IDS)
print(f"  emb={emb.shape}, ids={ids.shape}")

types = ids["type"].astype(str).values
unique_types = sorted(set(types))
type_to_idx = {t: i for i, t in enumerate(unique_types)}
type_oh = np.eye(len(unique_types), dtype=np.float32)[np.array([type_to_idx[t] for t in types])]
print(f"  types={unique_types}")

X = np.concatenate([emb[:].astype(np.float32), type_oh], axis=1)
del emb
print(f"  X={X.shape}, took {time.time()-t0:.1f}s")

# 用 step2 已训练好的 booster
print("\n加载 XGBoost 模型...")
booster = xgb.Booster()
booster.load_model(str(MODEL))
print(f"  best_iter={booster.best_iteration}")

# 打分
print("\n预测...")
dmat = xgb.DMatrix(X)
preds = booster.predict(dmat, iteration_range=(0, booster.best_iteration + 1))
print(f"  预测耗时 {time.time()-t0:.1f}s, preds shape={preds.shape}")

ids["pred_brightness_xgb"] = preds
ids["pred_relative_xgb"] = 10 ** preds  # brightness 是 log10,反 log 后是相对量级

# 整体 top 50
print("\n=== 141K 全库 XGBoost 预测 Top 50(整体) ===")
top50 = ids.nlargest(50, "pred_brightness_xgb")[["type", "mutations", "brightness", "pred_brightness_xgb", "pred_relative_xgb"]]
print(top50.to_string())

# 按 type 分组看 top
print("\n=== 各类型 Top 10 (按 XGBoost 预测) ===")
per_type_top = {}
for t in unique_types:
    sub = ids[ids["type"] == t].nlargest(10, "pred_brightness_xgb")
    per_type_top[t] = sub[["mutations", "brightness", "pred_brightness_xgb", "pred_relative_xgb"]]
    print(f"\n[{t}]")
    print(sub.to_string())

# 保存
ids.to_csv(ROOT / "step3a_141k_scored.csv", index=False)
summary = {
    "n_total": int(len(ids)),
    "best_overall": {
        "type": top50.iloc[0]["type"],
        "mutations": top50.iloc[0]["mutations"],
        "pred_brightness": float(top50.iloc[0]["pred_brightness_xgb"]),
        "pred_relative": float(top50.iloc[0]["pred_relative_xgb"]),
    },
    "per_type_best": {
        t: {
            "mutations": per_type_top[t].iloc[0]["mutations"],
            "pred_brightness": float(per_type_top[t].iloc[0]["pred_brightness_xgb"]),
            "pred_relative": float(per_type_top[t].iloc[0]["pred_relative_xgb"]),
        } for t in unique_types
    },
    "top50_count_by_type": top50["type"].value_counts().to_dict(),
}
with open(ROOT / "step3a_summary.json", "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)
print(f"\n=== Step 3a DONE ({time.time()-t0:.1f}s) ===")