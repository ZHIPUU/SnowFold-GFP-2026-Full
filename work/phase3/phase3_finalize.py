"""
Phase 3 Final: 终选 6 条 + 验证 + 生成 submission
==================================================
策略:
  - 5 条来自 phase2_scored.csv 的 ESM + 加性综合最优 (跨母体多样性)
  - 1 条 avGFP + sfGFP 经典超稳突变 (热稳定性保险)
  - 全部通过: 长度 220-250, M开头, exclusion list 比对
"""
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import esm

WORK = Path(r"D:\生信\2026Protein Design\work")
PHASE1 = WORK / "phase1"
PHASE2 = WORK / "phase2"
PHASE3 = WORK / "phase3"

# 加载 phase2 评分结果
cands = pd.read_csv(PHASE2 / "phase2_scored.csv")
print(f"Loaded {len(cands)} scored candidates")

# 加载 exclusion
excl_df = pd.read_csv(r"D:\生信\2026Protein Design\Exclusion_List.csv")
excl_seqs = set(excl_df["Sequence"].values)
print(f"Exclusion list: {len(excl_seqs)}")

# 加载 WT
with open(PHASE1 / "phase1_cache_v2.pkl", "rb") as f:
    cache = pickle.load(f)
wt_seqs = cache["wt_seqs"]
wt_brightness = cache["wt_brightness"]

# ---------- 选 5 条 (跨母体, ESM + 加性最优) ----------
# 排序每个 type
top_per_type = {}
for t in ["avGFP", "amacGFP", "cgreGFP", "ppluGFP"]:
    sub = cands[cands["type"] == t].sort_values("combined", ascending=False)
    top_per_type[t] = sub.head(3)  # 取每 type top 3

# 选 5 条:
# 1. avGFP best (2-mut 简洁)
# 2. avGFP 3-mut (备选)
# 3. cgreGFP best 4-mut (高加性预测)
# 4. ppluGFP best 4-mut (pplu最高分)
# 5. amacGFP best 2-mut (amac多样性)
picks_data = []
picks_data.append(top_per_type["avGFP"].iloc[0])  # avGFP best
picks_data.append(top_per_type["avGFP"].iloc[1])  # avGFP 2nd
picks_data.append(top_per_type["cgreGFP"].iloc[1])  # cgreGFP 4-mut best
picks_data.append(top_per_type["ppluGFP"].iloc[0])  # ppluGFP best
picks_data.append(top_per_type["amacGFP"].iloc[0])  # amacGFP best

# ---------- 构造第 6 条: avGFP + sfGFP 经典超稳突变 ----------
# avGFP 上的 sfGFP 经典位点都是有效的 (S65, S72, V163, T203, S202, A206, I171)
# 这些突变在加性模型里 NOT IN MODEL, 但物理上已知增强稳定性
# 我们额外加 1-2 个加性 top 突变 (I152S, E171Q)

# 加载 ESM2 用于最终验证
device = "cuda"
print("\nLoading ESM2-150M...")
model, alphabet = esm.pretrained.esm2_t30_150M_UR50D()
batch_converter = alphabet.get_batch_converter()
model = model.to(device).eval()

def compute_pll(seqs, batch_size=8):
    n = len(seqs)
    pll_per_res = np.zeros(n)
    pad_idx = alphabet.padding_idx
    eos_idx = alphabet.eos_idx
    bos_idx = alphabet.cls_idx
    unk_idx = alphabet.unk_idx
    for batch_start in range(0, n, batch_size):
        batch_seqs = seqs[batch_start:batch_start+batch_size]
        data = [(f"s_{i}", s) for i, s in enumerate(batch_seqs)]
        _, _, batch_tokens = batch_converter(data)
        batch_tokens = batch_tokens.to(device)
        with torch.no_grad():
            logits = model(batch_tokens)["logits"]
            log_probs = torch.log_softmax(logits, dim=-1)
            ll = log_probs.gather(2, batch_tokens.unsqueeze(-1)).squeeze(-1)
            mask = (
                (batch_tokens != pad_idx) &
                (batch_tokens != eos_idx) &
                (batch_tokens != bos_idx) &
                (batch_tokens != unk_idx)
            )
            seq_ll = (ll * mask).sum(dim=1).cpu().numpy()
            n_valid = mask.sum(dim=1).cpu().numpy()
        for i in range(len(batch_seqs)):
            pll_per_res[batch_start + i] = seq_ll[i] / max(n_valid[i], 1)
    return pll_per_res

def apply_mutations(wt_seq, muts):
    """muts: list of (pos, wt_aa, new_aa)"""
    seq = list(wt_seq)
    for pos, wt_aa, new_aa in muts:
        idx = pos - 1
        if idx < 0 or idx >= len(seq):
            return None
        if seq[idx] != wt_aa:
            return None
        seq[idx] = new_aa
    return "".join(seq)

# sfGFP 经典 + 加性 top 在 avGFP 上
sf_classical_avGFP = [
    (65, "S", "T"),    # S65T 发色团成熟
    (72, "S", "A"),    # S72A loop
    (163, "V", "A"),   # V163A β桶
    (203, "T", "I"),   # T203I β桶核心
    (202, "S", "D"),   # S202D 疏水核心
    (206, "A", "V"),   # A206V monomer
]

# 验证这些位置在 avGFP 上确实是这些 AA
avGFP = wt_seqs["avGFP"]
for pos, wt_aa, new_aa in sf_classical_avGFP:
    actual = avGFP[pos-1]
    print(f"  avGFP pos {pos}: expected {wt_aa}, actual {actual}, match: {actual == wt_aa}")

# 构造 6 号候选: avGFP + sfGFP核心 + I152S (加性 top 1)
mut6 = sf_classical_avGFP + [(152, "I", "S")]  # 7 突变
seq6 = apply_mutations(avGFP, mut6)
print(f"\nSeq 6 ({avGFP[:10]}... mutated): len={len(seq6)}, starts_with_M={seq6.startswith('M')}")
print(f"  in_exclusion: {seq6 in excl_seqs}")
print(f"  mutations: {'/'.join(f'{wt_aa}{pos}{new_aa}' for pos, wt_aa, new_aa in mut6)}")

# ---------- 验证 6 条 ----------
print("\n=== Validating final 6 ===")
final = []

# 5 条 from phase2
for i, row in enumerate(picks_data, 1):
    seq = row["seq"]
    name = f"S{i}_{row['type']}_{row['n_mut']}mut"
    final.append({
        "Seq_ID": i,
        "Name": name,
        "Type": row["type"],
        "n_mut": row["n_mut"],
        "mutations": row["mut_str"],
        "pred_brightness": row["pred_brightness"],
        "pred_relative": row["pred_relative"],
        "esm_pll_per_res": row["esm_pll_per_res"],
        "esm_delta_pll": row["esm_delta_pll"],
        "combined": row["combined"],
        "seq": seq,
        "len": len(seq),
        "in_exclusion": seq in excl_seqs,
        "starts_M": seq.startswith("M"),
    })

# 6 号 - sfGFP-classical
final.append({
    "Seq_ID": 6,
    "Name": "S6_avGFP_sfGFP_superfold",
    "Type": "avGFP",
    "n_mut": len(mut6),
    "mutations": "/".join(f"{wt_aa}{pos}{new_aa}" for pos, wt_aa, new_aa in mut6),
    "pred_brightness": None,  # sfGFP 突变不在模型里
    "pred_relative": None,
    "esm_pll_per_res": None,  # 待计算
    "esm_delta_pll": None,
    "combined": None,
    "seq": seq6,
    "len": len(seq6),
    "in_exclusion": seq6 in excl_seqs,
    "starts_M": seq6.startswith("M"),
})

# ESM 评估 Seq 6
seqs_for_esm = [final[5]["seq"], wt_seqs["avGFP"]]
plls = compute_pll(seqs_for_esm)
final[5]["esm_pll_per_res"] = plls[0]
final[5]["esm_delta_pll"] = plls[0] - plls[1]

# 输出
print(f"\n{'#':<3} {'Type':<10} {'nMut':<5} {'Len':<4} {'M':<3} {'Excl':<5} {'pred_b':<8} {'rel×':<7} {'PLL/res':<10} {'ΔPLL':<10} {'Description':<50}")
print("-" * 130)
for f in final:
    pb = f"{f['pred_brightness']:.3f}" if f['pred_brightness'] is not None else "N/A"
    pr = f"{f['pred_relative']:.2f}" if f['pred_relative'] is not None else "N/A"
    pl = f"{f['esm_pll_per_res']:.4f}" if f['esm_pll_per_res'] is not None else "N/A"
    dl = f"{f['esm_delta_pll']:+.4f}" if f['esm_delta_pll'] is not None else "N/A"
    print(f"{f['Seq_ID']:<3} {f['Type']:<10} {f['n_mut']:<5} {f['len']:<4} {str(f['starts_M'])[0]:<3} {str(f['in_exclusion'])[0]:<5} {pb:<8} {pr:<7} {pl:<10} {dl:<10} {f['mutations'][:50]:<50}")

# 全部检查
all_pass = all(220 <= f["len"] <= 250 and f["starts_M"] and not f["in_exclusion"] for f in final)
print(f"\nAll candidates pass constraints: {all_pass}")

# 保存
df_final = pd.DataFrame(final)
df_final.to_csv(PHASE3 / "final_6_candidates.csv", index=False)
print(f"\n[OK] Saved to {PHASE3 / 'final_6_candidates.csv'}")

# ---------- 生成 submission csv ----------
sub_df = pd.DataFrame({
    "Team_Name": ["YourTeamName"] * 6,
    "Seq_ID": [f["Seq_ID"] for f in final],
    "Sequence": [f["seq"] for f in final],
})
sub_df.to_csv(WORK / "submission_yourteamname.csv", index=False)
print(f"[OK] Submission csv saved to {WORK / 'submission_yourteamname.csv'}")

# ---------- 输出每条序列的最终信息 ----------
print("\n" + "=" * 80)
print("FINAL 6 CANDIDATES - DETAILED")
print("=" * 80)
for f in final:
    print(f"\nSeq {f['Seq_ID']}: {f['Name']}")
    print(f"  Type: {f['Type']} ({len(wt_seqs[f['Type']])} aa WT)")
    print(f"  Mutations ({f['n_mut']}): {f['mutations']}")
    print(f"  Length: {f['len']} aa")
    print(f"  Starts with M: {f['starts_M']}")
    print(f"  In exclusion list: {f['in_exclusion']}")
    if f['pred_brightness'] is not None:
        print(f"  Additive model pred brightness: {f['pred_brightness']:.3f} (rel {f['pred_relative']:.2f}× WT)")
    if f['esm_pll_per_res'] is not None:
        wt_pll = compute_pll([wt_seqs[f['Type']]])[0]
        print(f"  ESM2 PLL/res: {f['esm_pll_per_res']:.4f} (vs WT {wt_pll:.4f}, Δ={f['esm_delta_pll']:+.4f})")
    print(f"  Sequence:")
    for i in range(0, len(f['seq']), 60):
        print(f"    {f['seq'][i:i+60]}")

print("\n=== Phase 3+4 DONE ===")