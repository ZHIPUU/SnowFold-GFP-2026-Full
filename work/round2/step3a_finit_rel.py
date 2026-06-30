"""Step 3a-extended: 计算每个 type 的 max Finit/Finit_WT,看数据天花板。"""
import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design\work\round2")
ids = pd.read_csv(ROOT / "step3a_141k_scored.csv")

WT = {"avGFP": 3.72, "amacGFP": 3.97, "cgreGFP": 4.50, "ppluGFP": 4.23}

ids["finit_rel_xgb"] = ids.apply(lambda r: 10 ** (r["pred_brightness_xgb"] - WT[r["type"]]), axis=1)

print("=== 各 type 的 Finit/Finit_WT 最大值(基于 XGBoost 预测) ===")
for t in ["avGFP", "amacGFP", "cgreGFP", "ppluGFP"]:
    sub = ids[ids["type"] == t]
    top = sub.nlargest(5, "finit_rel_xgb")[["mutations", "brightness", "pred_brightness_xgb", "finit_rel_xgb"]]
    print(f"\n[{t}] WT baseline log10={WT[t]}, Finit_rel 最大:")
    print(top.to_string())

print("\n=== 整体 Finit/Finit_WT 分布 ===")
ids["log10_finit_rel"] = np.log10(ids["finit_rel_xgb"])
print(ids.groupby("type")["log10_finit_rel"].describe().round(3).to_string())

# max over all types
max_rel = ids["finit_rel_xgb"].max()
print(f"\n所有 type 中 Finit_rel 最大: {max_rel:.2f}×")
print(f"对应 brightness 预测: {ids.loc[ids['finit_rel_xgb'].idxmax(), 'pred_brightness_xgb']:.3f}")
print(f"对应 type: {ids.loc[ids['finit_rel_xgb'].idxmax(), 'type']}, mutations: {ids.loc[ids['finit_rel_xgb'].idxmax(), 'mutations']}")

# 实际测量的 max
ids["finit_rel_actual"] = ids.apply(lambda r: 10 ** (r["brightness"] - WT[r["type"]]), axis=1)
max_actual = ids["finit_rel_actual"].max()
print(f"\n实际测量的 Finit_rel 最大: {max_actual:.2f}×")
print(f"对应 brightness: {ids.loc[ids['finit_rel_actual'].idxmax(), 'brightness']:.3f}")
print(f"对应 type: {ids.loc[ids['finit_rel_actual'].idxmax(), 'type']}, mutations: {ids.loc[ids['finit_rel_actual'].idxmax(), 'mutations']}")