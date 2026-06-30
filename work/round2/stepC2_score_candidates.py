"""Step C2 (修正): 用 epistasis 模型重打分 11 候选,选 Top-6 更新 submission。

正确做法: 从 WT vs candidate 序列对比提取突变,再算 epistasis 特征。"""
import time
import re
import json
import numpy as np
import pandas as pd
import xgboost as xgb
from itertools import combinations
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design\work\round2")
MODEL = ROOT / "stepC_xgboost_epistasis.model"
CAND_FILE = ROOT / "candidates_round2_design.csv"
EMB_PATH = ROOT / "esm650m_embeddings.npy"
IDS_PATH = ROOT / "esm650m_ids.csv"
WT_FILE = Path(r"D:\生信\2026Protein Design\AAseqs of 5 GFP proteins_20260511.txt")
OUT_SUB = Path(r"D:\生信\2026Protein Design\submission_yourteamname.csv")

t0 = time.time()

# 加载 WT
wt = {}
current_name = None
with open(WT_FILE) as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            current_name = line[1:]
            wt[current_name] = ""
        elif line and not line.startswith("#"):
            wt[current_name] += line

WT = {"avGFP": 3.72, "amacGFP": 3.97, "cgreGFP": 4.50, "ppluGFP": 4.23}
# sfGFP 没有官方 baseline;近似用 avGFP WT + sfGFP 11 mut 经验增量 ~0.45 → 4.17
WT["sfGFP"] = 4.17

# 加载 epistasis 模型
print("加载 epistasis 模型...")
booster = xgb.Booster()
booster.load_model(str(MODEL))
print(f"  best_iter={booster.best_iteration}")

# 加载 141K 数据,估计 mutation effects
print("加载 141K 数据...")
emb_full = np.load(EMB_PATH, mmap_mode="r")
ids_full = pd.read_csv(IDS_PATH)

def parse_muts(s):
    if pd.isna(s) or not s: return []
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

ids_full["mut_list"] = ids_full["mutations"].apply(parse_muts)
ids_full["n_mut"] = ids_full["mut_list"].apply(len)

print("估计单点效应...")
single_effects_mean = {}
for parent in WT:
    sub = ids_full[ids_full["type"] == parent]
    single_subs = sub[sub["n_mut"] == 1]
    for _, row in single_subs.iterrows():
        frm, pos, to = row["mut_list"][0]
        key = (parent, frm, pos, to)
        delta = row["brightness"] - WT[parent]
        single_effects_mean.setdefault(key, []).append(delta)
single_effects_mean = {k: np.mean(v) for k, v in single_effects_mean.items()}
print(f"  {len(single_effects_mean)} unique single effects")

print("估计 pairwise epistasis...")
pair_epistasis_mean = {}
for parent in WT:
    sub = ids_full[ids_full["type"] == parent]
    pair_subs = sub[sub["n_mut"] == 2]
    for _, row in pair_subs.iterrows():
        muts = row["mut_list"]
        additive = 0.0
        valid = True
        for frm, pos, to in muts:
            if (parent, frm, pos, to) in single_effects_mean:
                additive += single_effects_mean[(parent, frm, pos, to)]
            else:
                valid = False
                break
        if not valid: continue
        epsilon = (row["brightness"] - WT[parent]) - additive
        key2 = (parent, frozenset([(frm, pos, to) for frm, pos, to in muts]))
        pair_epistasis_mean.setdefault(key2, []).append(epsilon)
pair_epistasis_mean = {k: np.mean(v) for k, v in pair_epistasis_mean.items()}
print(f"  {len(pair_epistasis_mean)} unique pair epistasis")

# 嵌入候选
print("\n嵌入 11 候选序列...")
import torch, esm
ckpt = torch.load(r"C:\Users\A\.cache\torch\hub\checkpoints\esm2_t33_650M_UR50D.pt.75c71e41769c4391ba0186bc8c92d0f7.partial",
                  map_location="cuda", weights_only=False)
args = ckpt["cfg"]["model"]
alphabet = esm.Alphabet.from_architecture("ESM-1b")
model = esm.ESM2(
    num_layers=args.encoder_layers,
    embed_dim=args.encoder_embed_dim,
    attention_heads=args.encoder_attention_heads,
    alphabet=alphabet,
).cuda()
model.load_state_dict(ckpt["model"], strict=False)
model.eval()

cand = pd.read_csv(CAND_FILE)
BATCH = 32
batch_converter = alphabet.get_batch_converter()
emb_cand = np.zeros((len(cand), 1280), dtype=np.float32)
for i in range(0, len(cand), BATCH):
    seqs = cand["seq"].iloc[i:i+BATCH].tolist()
    data = [(f"p{j}", s[:1022]) for j, s in enumerate(seqs)]
    _, _, bt = batch_converter(data)
    bt = bt.cuda()
    with torch.no_grad():
        out = model(bt, repr_layers=[33], return_contacts=False)
    emb_cand[i:i+len(seqs)] = out["representations"][33][:, 1:-1].mean(dim=1).cpu().numpy()
print(f"  嵌入 shape: {emb_cand.shape}, took {time.time()-t0:.1f}s")
del model
torch.cuda.empty_cache()

# 计算候选的 epistasis 特征(从 scaffold WT 对比)
print("\n从 scaffold WT 对比提取候选突变...")
unique_types = ['amacGFP', 'avGFP', 'cgreGFP', 'ppluGFP']
type_to_idx = {t: i for i, t in enumerate(unique_types)}
# 处理未知 scaffold (如 sfGFP) -> 用零向量
scaffold_idx = cand["scaffold"].map(type_to_idx)
scaffold_idx = scaffold_idx.fillna(-1).astype(int).values
type_oh = np.zeros((len(cand), len(unique_types)), dtype=np.float32)
valid_mask = scaffold_idx >= 0
type_oh[valid_mask] = np.eye(len(unique_types), dtype=np.float32)[scaffold_idx[valid_mask]]
print(f"  type_oh shape: {type_oh.shape}, valid scaffolds: {valid_mask.sum()}/{len(cand)}")

single_sum = np.zeros(len(cand), dtype=np.float32)
single_mean = np.zeros(len(cand), dtype=np.float32)
single_max = np.zeros(len(cand), dtype=np.float32)
single_min = np.zeros(len(cand), dtype=np.float32)
pair_eps_sum = np.zeros(len(cand), dtype=np.float32)
pair_eps_count = np.zeros(len(cand), dtype=np.float32)
unseen_pair_count = np.zeros(len(cand), dtype=np.float32)
n_mut_arr = np.zeros(len(cand), dtype=np.float32)

mut_lists = []
for i, row in cand.iterrows():
    parent = row["scaffold"]
    wt_seq = wt[parent]
    cand_seq = row["seq"]
    # 对比提取突变
    muts = []
    # 处理长度差异(cgreGFP/ppluGFP 长度不同)
    min_len = min(len(wt_seq), len(cand_seq))
    for pos in range(min_len):
        if wt_seq[pos] != cand_seq[pos]:
            muts.append((wt_seq[pos], pos+1, cand_seq[pos]))
    mut_lists.append(muts)
    n_mut_arr[i] = len(muts)
    if not muts:
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

# 拼接特征
extra = np.column_stack([
    type_oh,  # 4
    n_mut_arr.reshape(-1, 1),  # 1
    single_sum.reshape(-1, 1),  # 1
    single_mean.reshape(-1, 1),  # 1
    single_max.reshape(-1, 1),  # 1
    single_min.reshape(-1, 1),  # 1
    pair_eps_sum.reshape(-1, 1),  # 1
    pair_eps_count.reshape(-1, 1),  # 1
    unseen_pair_count.reshape(-1, 1),  # 1
])
X_cand = np.concatenate([emb_cand, extra], axis=1)
print(f"  X_cand shape: {X_cand.shape}")

# 预测
print("\n预测...")
dmat = xgb.DMatrix(X_cand)
preds = booster.predict(dmat, iteration_range=(0, booster.best_iteration + 1))
cand["pred_brightness"] = preds
cand["finit_rel"] = [10 ** (preds[i] - WT[cand["scaffold"].iloc[i]]) for i in range(len(cand))]
cand["n_mut_calc"] = [len(m) for m in mut_lists]

print("\n=== 候选按 epistasis 模型排序 ===")
print(cand[["id","name","scaffold","n_mut_calc","pred_brightness","finit_rel"]].sort_values("finit_rel", ascending=False).to_string(index=False))

cand.to_csv(ROOT / "stepC2_candidates_scored.csv", index=False)

# 选 Top-6(平衡多样性 + 高 score)
sorted_cand = cand.sort_values("finit_rel", ascending=False)
selected = []
selected_scaffolds = []
# 策略: 优先 score 高,然后 scaffold 多样
for _, row in sorted_cand.iterrows():
    if len(selected) >= 6:
        break
    sc = row["scaffold"]
    # 限制每 scaffold 最多 2 条
    n_sc = selected_scaffolds.count(sc)
    if n_sc >= 2:
        continue
    selected.append(int(row["id"]))
    selected_scaffolds.append(sc)

print(f"\n=== Top-6 selected: {selected} (scaffolds: {selected_scaffolds}) ===")

# 生成 submission
rows = []
for new_id, cid in enumerate(selected, 1):
    row = cand[cand["id"]==cid].iloc[0]
    rows.append({"Team_Name": "YourTeamName", "Seq_ID": new_id, "Sequence": row["seq"]})
df_sub = pd.DataFrame(rows)
df_sub.to_csv(OUT_SUB, index=False)
print(f"\n✅ Submission 已更新到 {OUT_SUB}")

# 验证
print("\n=== 最终验证 ===")
for new_id, cid in enumerate(selected, 1):
    row = cand[cand["id"]==cid].iloc[0]
    s = row["seq"]
    chrom = "TYG" in s or "SYG" in s or "GYG" in s
    print(f"  Seq {new_id}: {row['name']:35s} len={len(s)} finit_rel={row['finit_rel']:.3f}× chromo={'✓' if chrom else '✗'}")