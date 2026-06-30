"""Step 3b: ESM2-650M 重嵌入 top 10K 候选(按老 Ridge 模型筛选) + XGBoost 重新打分。"""
import time
import json
import numpy as np
import pandas as pd
import torch
import xgboost as xgb
import esm
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design\work\round2")
TOP_CAND = Path(r"D:\生信\2026Protein Design\work\phase1\top_candidates.csv")
ESM_MODEL_PATH = Path(r"C:\Users\A\.cache\torch\hub\checkpoints\esm2_t33_650M_UR50D.pt.75c71e41769c4391ba0186bc8c92d0f7.partial")

t0 = time.time()
print("加载 top_candidates.csv ...")
cand = pd.read_csv(TOP_CAND)
print(f"  total {len(cand)} rows, cols={list(cand.columns)}")

# 按旧 Ridge 模型 pred_relative 排序,取 top 10K
cand = cand.sort_values("pred_relative", ascending=False).head(10000).reset_index(drop=True)
print(f"  top 10K by old Ridge, pred_relative range [{cand['pred_relative'].min():.2f}, {cand['pred_relative'].max():.2f}]")

# 加载 ESM2-650M
print("\n加载 ESM2-650M ...")
ckpt = torch.load(str(ESM_MODEL_PATH), map_location="cuda", weights_only=False)
args = ckpt["cfg"]["model"]
alphabet = esm.Alphabet.from_architecture("ESM-1b")
model = esm.ESM2(
    num_layers=args.encoder_layers,
    embed_dim=args.encoder_embed_dim,
    attention_heads=args.encoder_attention_heads,
    alphabet=alphabet,
).cuda()
missing, unexpected = model.load_state_dict(ckpt["model"], strict=False)
model.eval()
print(f"  加载完毕, missing={len(missing)}, unexpected={len(unexpected)}")

# 批量嵌入
BATCH = 128
batch_converter = alphabet.get_batch_converter()
embeddings = np.zeros((len(cand), 1280), dtype=np.float32)

t1 = time.time()
for i in range(0, len(cand), BATCH):
    seqs = cand["seq"].iloc[i:i+BATCH].tolist()
    data = [(f"p{j}", s[:1022]) for j, s in enumerate(seqs)]
    _, _, batch_tokens = batch_converter(data)
    batch_tokens = batch_tokens.cuda()
    with torch.no_grad():
        out = model(batch_tokens, repr_layers=[33], return_contacts=False)
    emb = out["representations"][33][:, 1:-1].mean(dim=1).cpu().numpy()
    embeddings[i:i+len(seqs)] = emb
    if (i // BATCH) % 5 == 0:
        elapsed = time.time() - t1
        rate = (i + len(seqs)) / max(elapsed, 0.01)
        eta = (len(cand) - i - len(seqs)) / max(rate, 0.01)
        print(f"  [{time.strftime('%H:%M:%S')}] {i+len(seqs)}/{len(cand)} ({100*(i+len(seqs))/len(cand):.1f}%) "
              f"rate={rate:.1f} seq/s, ETA={eta:.0f}s, GPU mem={torch.cuda.memory_allocated()/1e9:.2f} GB")

print(f"\n嵌入完成, 耗时 {time.time()-t1:.1f}s, embeddings shape={embeddings.shape}")
np.save(ROOT / "top10k_esm650m_embeddings.npy", embeddings)
cand.to_csv(ROOT / "top10k_candidates.csv", index=False)

# 用 XGBoost GPU 打分
print("\n加载 XGBoost 模型 ...")
booster = xgb.Booster()
booster.load_model(str(ROOT / "step2_xgboost_gpu.model"))

# 构建特征
unique_types = ['amacGFP', 'avGFP', 'cgreGFP', 'ppluGFP']
type_to_idx = {t: i for i, t in enumerate(unique_types)}
type_oh = np.eye(len(unique_types), dtype=np.float32)[np.array([type_to_idx[t] for t in cand["type"]])]
X = np.concatenate([embeddings, type_oh], axis=1)
print(f"  X={X.shape}")

dmat = xgb.DMatrix(X)
preds = booster.predict(dmat, iteration_range=(0, booster.best_iteration + 1))
print(f"  XGBoost 预测完成, shape={preds.shape}")

WT = {"avGFP": 3.72, "amacGFP": 3.97, "cgreGFP": 4.50, "ppluGFP": 4.23}
cand["pred_brightness_xgb"] = preds
cand["finit_rel_xgb"] = [10 ** (preds[j] - WT[cand["type"].iloc[j]]) for j in range(len(cand))]
cand.to_csv(ROOT / "top10k_scored.csv", index=False)

# 输出整体 + 各类型 top
print("\n=== Top 50 (按 XGBoost Finit_rel) ===")
top50 = cand.nlargest(50, "finit_rel_xgb")[["type", "mut_str", "pred_brightness_xgb", "finit_rel_xgb"]]
print(top50.to_string())

per_type = {}
for t in unique_types:
    sub = cand[cand["type"] == t]
    if len(sub) == 0:
        continue
    per_type[t] = {
        "n": int(len(sub)),
        "best_mut": sub.nlargest(1, "finit_rel_xgb").iloc[0]["mut_str"],
        "best_finit_rel": float(sub["finit_rel_xgb"].max()),
        "top10": sub.nlargest(10, "finit_rel_xgb")[["mut_str", "finit_rel_xgb"]].to_dict(orient="records"),
    }

with open(ROOT / "step3b_summary.json", "w", encoding="utf-8") as f:
    json.dump(per_type, f, indent=2, ensure_ascii=False)

print("\n=== 各 type 最佳 Finit/Finit_WT ===")
for t in unique_types:
    if t in per_type:
        print(f"  [{t}] best_mut={per_type[t]['best_mut']}, finit_rel={per_type[t]['best_finit_rel']:.3f}×")

print(f"\n=== Step 3b DONE (总耗时 {time.time()-t0:.1f}s) ===")