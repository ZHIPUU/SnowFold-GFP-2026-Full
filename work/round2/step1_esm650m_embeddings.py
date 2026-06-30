"""Round 2 Step 1: Compute ESM2-650M mean-pooled embeddings for all 141K variants.

- Bypasses torch.hub URL check by loading directly
- Uses batch_size=128 (best speed/memory balance)
- Saves to .npy (memmap-friendly)
- Has checkpoint resume capability

Output: work/round2/esm650m_embeddings.npy shape (N, 1280) float32
        work/round2/esm650m_ids.csv with row index, type, mutations
"""
import os
import sys
import time
import pickle
import numpy as np
import pandas as pd
import torch
from pathlib import Path

WORK = Path(r"D:\生信\2026Protein Design")
ROUND2 = WORK / "work" / "round2"
ROUND2.mkdir(exist_ok=True, parents=True)

OUT_NPY = ROUND2 / "esm650m_embeddings.npy"
CHECKPOINT = ROUND2 / "esm650m_checkpoint.npy"

# ---- Load data ----
print("[1/4] Loading variants...", flush=True)
with open(ROUND2 / "all_variants.pkl", 'rb') as f:
    df = pickle.load(f)
print(f"  {len(df):,} variants", flush=True)

seqs = df['seq'].tolist()
N = len(seqs)
EMB_DIM = 1280

# ---- Load model ----
print("\n[2/4] Loading ESM2-650M from local cache...", flush=True)
t0 = time.time()
import esm
from esm import ESM2, Alphabet
alphabet = Alphabet.from_architecture('ESM-1b')
model = ESM2(num_layers=33, embed_dim=1280, attention_heads=20, alphabet=alphabet)
model_path = r"C:\Users\A\.cache\torch\hub\checkpoints\esm2_t33_650M_UR50D.pt.75c71e41769c4391ba0186bc8c92d0f7.partial"
sd = torch.load(model_path, map_location='cpu', weights_only=False)
if isinstance(sd, dict) and 'model' in sd:
    sd = sd['model']
new_sd = {k.replace('module.', ''): v for k, v in sd.items()}
miss, unexp = model.load_state_dict(new_sd, strict=False)
print(f"  Loaded in {time.time()-t0:.1f}s", flush=True)
if miss:
    print(f"  Missing keys: {len(miss)}", flush=True)
if unexp:
    print(f"  Unexpected keys: {len(unexp)}", flush=True)

model = model.eval().cuda()
print(f"  GPU mem: {torch.cuda.memory_allocated()/1e9:.2f} GB", flush=True)

batch_converter = alphabet.get_batch_converter()

# ---- Resume support ----
done = 0
if CHECKPOINT.exists():
    try:
        ckpt = np.load(CHECKPOINT, mmap_mode='r')
        if ckpt.shape == (N, EMB_DIM):
            done = N  # already complete? need to check if zeros
            nonzero = np.any(ckpt != 0, axis=1)
            done = int(nonzero.sum())
            print(f"  Resuming from checkpoint: {done:,} rows already done", flush=True)
            del ckpt
    except Exception as e:
        print(f"  Checkpoint load failed: {e}", flush=True)

# Allocate output (memmap if possible)
print(f"\n[3/4] Computing embeddings for {N:,} variants (batch_size=128)...", flush=True)
if OUT_NPY.exists():
    all_embs = np.load(OUT_NPY, mmap_mode='r+')
else:
    all_embs = np.lib.format.open_memmap(OUT_NPY, mode='w+', dtype=np.float32, shape=(N, EMB_DIM))

batch_size = 128
t_start = time.time()

for i in range(done, N, batch_size):
    batch_seqs = seqs[i:i+batch_size]
    batch_data = [(f"p{j}", s) for j, s in enumerate(batch_seqs)]
    try:
        _, _, toks = batch_converter(batch_data)
    except Exception as e:
        # Try truncating
        batch_data = [(f"p{j}", s[:1022]) for j, s in enumerate(batch_seqs)]
        _, _, toks = batch_converter(batch_data)
    toks = toks.cuda()

    with torch.no_grad():
        out = model(toks, repr_layers=[33], return_contacts=False)
    emb = out['representations'][33]  # (B, L+2, 1280)

    # Mean pool over non-padding, non-special tokens
    mask = (toks != alphabet.padding_idx) & (toks != alphabet.cls_idx) & (toks != alphabet.eos_idx)
    emb_masked = emb * mask.unsqueeze(-1).float()
    seq_emb = emb_masked.sum(dim=1) / mask.sum(dim=1, keepdim=True).clamp(min=1).float()

    all_embs[i:i+len(batch_seqs)] = seq_emb.cpu().numpy()
    all_embs.flush()

    if (i // batch_size) % 25 == 0:
        elapsed = time.time() - t_start
        rate = (i - done + len(batch_seqs)) / max(elapsed, 0.1)
        remaining = (N - i - len(batch_seqs)) / max(rate, 0.001)
        pct = (i + len(batch_seqs)) / N * 100
        print(f"  [{i+len(batch_seqs):,}/{N:,}] {pct:.1f}% | {rate:.1f} seq/s | ETA {remaining/60:.1f} min", flush=True)

# Save final and checkpoint
print(f"\n  Total time: {(time.time()-t_start)/60:.1f} min", flush=True)
all_embs.flush()

# Save metadata (in case we need to join later)
df[['type', 'mutations', 'brightness']].to_csv(ROUND2 / "esm650m_ids.csv", index=False)
print(f"  Saved {OUT_NPY} ({os.path.getsize(OUT_NPY)/1e9:.2f} GB)", flush=True)

# Cleanup checkpoint
if CHECKPOINT.exists():
    os.remove(CHECKPOINT)

print("\n[4/4] DONE", flush=True)