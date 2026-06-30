"""
Phase 2: ESM2 候选打分 (高效版)
==================================================
方法: 一次 forward, 取每个位置真实 AA 的 log-prob (pseudo-log-likelihood)
对每个序列返回总 log-likelihood 和 per-residue 平均
"""
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import esm

WORK = Path(r"D:\生信\2026Protein Design\work")
PHASE1 = WORK / "phase1"
PHASE2 = WORK / "phase2"
PHASE2.mkdir(parents=True, exist_ok=True)

device = "cuda"
print("Loading ESM2-150M...")
t0 = time.time()
model, alphabet = esm.pretrained.esm2_t30_150M_UR50D()
batch_converter = alphabet.get_batch_converter()
model = model.to(device).eval()
print(f"  Load: {time.time()-t0:.1f}s, GPU mem: {torch.cuda.memory_allocated()/1024**3:.2f} GB")

# ---------- 加载候选 ----------
cands = pd.read_csv(PHASE1 / "top50_candidates.csv")
print(f"Loaded {len(cands)} candidates")

# ---------- 计算 pseudo-log-likelihood ----------
def compute_pll(seqs, model, batch_converter, device, batch_size=16):
    """
    一次 forward, 评估每个位置的 log-prob
    返回每个序列的总 PLL 和 per-residue PLL
    """
    n = len(seqs)
    pll_total = np.zeros(n)
    pll_per_res = np.zeros(n)

    # 关键 token ids
    pad_idx = alphabet.padding_idx
    eos_idx = alphabet.eos_idx
    bos_idx = alphabet.cls_idx  # BOS
    unk_idx = alphabet.unk_idx

    # 逐 batch 跑
    for batch_start in range(0, n, batch_size):
        batch_seqs = seqs[batch_start:batch_start+batch_size]
        data = [(f"s_{i}", s) for i, s in enumerate(batch_seqs)]
        _, _, batch_tokens = batch_converter(data)
        batch_tokens = batch_tokens.to(device)
        B, L = batch_tokens.shape

        with torch.no_grad():
            logits = model(batch_tokens)["logits"]  # [B, L, vocab]
            log_probs = torch.log_softmax(logits, dim=-1)

            # 取真实 token 的 log-prob
            ll = log_probs.gather(2, batch_tokens.unsqueeze(-1)).squeeze(-1)  # [B, L]

            # 有效 mask: 排除 BOS, EOS, padding, 未知
            mask = (
                (batch_tokens != pad_idx) &
                (batch_tokens != eos_idx) &
                (batch_tokens != bos_idx) &
                (batch_tokens != unk_idx)
            )

            seq_ll = (ll * mask).sum(dim=1).cpu().numpy()  # [B]
            n_valid = mask.sum(dim=1).cpu().numpy()  # [B]

        for i in range(len(batch_seqs)):
            idx = batch_start + i
            pll_total[idx] = seq_ll[i]
            pll_per_res[idx] = seq_ll[i] / max(n_valid[i], 1)

        if (batch_start // batch_size) % 5 == 0:
            print(f"  batch {batch_start//batch_size + 1}/{(n+batch_size-1)//batch_size} done")

    return pll_total, pll_per_res


# 跑打分
print("\n=== Computing ESM2 PLL for all candidates ===")
t0 = time.time()
pll_total, pll_per_res = compute_pll(cands["seq"].tolist(), model, batch_converter, device, batch_size=16)
print(f"  Total time: {time.time()-t0:.1f}s")
print(f"  PLL range: {pll_total.min():.2f} to {pll_total.max():.2f}")
print(f"  PLL/res range: {pll_per_res.min():.4f} to {pll_per_res.max():.4f}")

cands["esm_pll"] = pll_total
cands["esm_pll_per_res"] = pll_per_res

# 也对 WT 序列打分作为 baseline
print("\n=== WT baselines ===")
wt_seqs = {
    "avGFP": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
    "amacGFP": "MSKGEELFTGIVPVLIELDGDVHGHKFSVRGEGEGDADYGKLEIKFICTTGKLPVPWPTLVTTLSYGILCFARYPEHMKMNDFFKSAMPEGYIQERTIFFQDDGKYKTRGEVKFEGDTLVNRIELKGMDFKEDGNILGHKLEYNFNSHNVYIMPDKANNGLKVNFKIRHNIEGGGVQLADHYQTNVPLGDGPVLIPINHYLSCQTAISKDRNETRDHMVFLEFFSACGHTHGMDELYK",
    "cgreGFP": "MTALTEGAKLFEKEIPYITELEGDVEGMKFIIKGEGTGDATTGTIKAKYICTTGDLPVPWATILSSLSYGVFCFAKYPRHIADFFKSTQPDGYSQDRIISFDNDGQYDVKAKVTYENGTLYNRVTVKGTGFKSNGNILGMRVLYHSPPHAVYILPDRKNGGMKIEYNKAFDVMGGGHQMARHAQFNKPLGAWEEDYPLYHHLTVWTSFGKDPDDDETDHLTIVEVIKAVDLETYR",
    "ppluGFP": "MPAMKIECRITGTLNGVEFELVGGGEGTPEQGRMTNKMKSTKGALTFSPYLLSHVMGYGFYHFGTYPSGYENPFLHAINNGGYTNTRIEKYEDGGVLHVSFSYRYEAGRVIGDFKVVGTGFPEDSVIFTDKIIRSNATVEHLHPMGDNVLVGSFARTFSLRDGGYYSFVVDSHMHFKSAIHPSILQNGGPMFAFRRVEELHSNTELGIVEYQHAFKTPIAFA",
}
wt_total, wt_per_res = compute_pll(list(wt_seqs.values()), model, batch_converter, device, batch_size=4)
wt_pll_dict = dict(zip(wt_seqs.keys(), wt_per_res))
for (name, _), s, sr in zip(wt_seqs.items(), wt_total, wt_per_res):
    print(f"  {name} WT: PLL={s:.2f}, PLL/res={sr:.4f}")

# 跟 WT 比
cands["wt_pll_per_res"] = cands["type"].map(wt_pll_dict)
cands["esm_delta_pll"] = cands["esm_pll_per_res"] - cands["wt_pll_per_res"]

# 联合排序
print("\n=== Top candidates by combined score ===")
print("(Combined = pred_brightness + esm_delta_pll * scale)")
for scale in [50, 100, 200, 500]:
    cands["combined"] = cands["pred_brightness"] + cands["esm_delta_pll"] * scale
    top = cands.sort_values("combined", ascending=False).head(5)
    print(f"\n  Scale={scale}:")
    for _, r in top.iterrows():
        print(f"    {r['type']:8s} {r['n_mut']}-mut: pred={r['pred_brightness']:.3f}, "
              f"esm_d={r['esm_delta_pll']:+.4f}, combined={r['combined']:.3f}  {r['mut_str']}")

# 保存
cands.to_csv(PHASE2 / "phase2_scored.csv", index=False)
print(f"\n[OK] Saved to {PHASE2 / 'phase2_scored.csv'}")

# Top 30 (scale=200)
cands["combined"] = cands["pred_brightness"] + cands["esm_delta_pll"] * 200
cands_sorted = cands.sort_values("combined", ascending=False)
print("\n=== Final Top 30 (scale=200) ===")
for i, (_, r) in enumerate(cands_sorted.head(30).iterrows()):
    print(f"  {i+1:2d}. {r['type']:8s} {r['n_mut']}-mut: "
          f"pred_b={r['pred_brightness']:.3f}, esm_d={r['esm_delta_pll']:+.4f}, "
          f"combined={r['combined']:.3f}  {r['mut_str']}")

print("\n=== Phase 2 DONE ===")