"""EVOLVEpro 评分 (模型已训练, 仅推理候选)"""
import numpy as np
import pandas as pd
import json
import time
from pathlib import Path
import xgboost as xgb
import torch
from transformers import AutoTokenizer, AutoModel

ROOT = Path(r"D:\生信\2026Protein Design")
R2 = ROOT / "work" / "round2"
R5 = ROOT / "work" / "round5"

# 加载已训练的 XGBoost 模型
booster = xgb.Booster()
booster.load_model(str(R5 / "evolvepro_xgb.model"))
print("XGBoost 加载完成")

# 加载训练数据 (only metadata 用于 type 映射 + WT 参考)
ids = pd.read_csv(R2 / "esm650m_ids.csv")
types = ids["type"].astype(str).values
unique_types = sorted(set(types))
type_to_idx = {t: i for i, t in enumerate(unique_types)}
y = ids["brightness"].astype(np.float32).values

# 加载所有候选
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

all_c = r4_hc + r4_msf + r4_mav + r5_lm + r5_le
seen = set(); dedup = []
for c in all_c:
    if c["seq"] not in seen:
        seen.add(c["seq"]); dedup.append(c)
print(f"候选 (去重): {len(dedup)}")

# 加载 ESM-2 650M
print("加载 ESM-2 650M (GPU)...")
tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t33_650M_UR50D")
model = AutoModel.from_pretrained("facebook/esm2_t33_650M_UR50D").cuda().eval()
print("OK\n")

scaffold_to_type = {
    "sfGFP": "avGFP", "avGFP": "avGFP", "amacGFP": "amacGFP",
    "cgreGFP": "cgreGFP", "ppluGFP": "ppluGFP",
    "mBaoJin": "avGFP",
    "sfGFP_MPNN": "avGFP", "avGFP_MPNN": "avGFP", "avGFP_LMPNN": "avGFP",
}

t0 = time.time()
for i, c in enumerate(dedup):
    if i % 10 == 0:
        print(f"  {i}/{len(dedup)} ({time.time()-t0:.1f}s)")
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

print(f"\n完成: {time.time()-t0:.1f}s")

# 排序
sorted_c = sorted(dedup, key=lambda x: -x["evolvepro_pred"])

print("\n" + "=" * 100)
print("Top-20 候选 (按 EVOLVEpro 预测亮度):")
print("=" * 100)
print(f"{'#':<3} {'name':<35} {'scaffold':<14} {'pLDDT':>5} {'mut':>3} {'EVOLVEpro':>9}")
print("-" * 95)
for i, c in enumerate(sorted_c[:20], 1):
    print(f"{i:<3} {c['name'][:35]:<35} {c['scaffold'][:14]:<14} "
          f"{c['plddt_mean']:>5.1f} {c['n_muts']:>3} {c['evolvepro_pred']:>9.3f}")

print(f"\nWT 训练集 brightness 参考 (log10):")
for t in unique_types:
    t_mask = types == t
    print(f"  {t}: mean={y[t_mask].mean():.3f}, max={y[t_mask].max():.3f}")
print(f"  整体: mean={y.mean():.3f}, max={y.max():.3f}")

with open(R5 / "evolvepro_scored.json", "w", encoding="utf-8") as f:
    json.dump(sorted_c, f, indent=2, ensure_ascii=False)
print(f"\n保存 evolvepro_scored.json")
