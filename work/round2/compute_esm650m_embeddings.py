"""Round 2: Compute ESM2-650M embeddings for all 140K variants in GFP_data.xlsx.

Strategy:
1. Load ESM2-650M model directly from local file (avoid URL check)
2. For each variant, compute mean-pooled embedding (sequence-level representation)
3. Save to .npy for downstream ML training

Output: work/round2/esm650m_embeddings.npy  shape: (N, 1280)
        work/round2/variant_ids.parquet (sequence + type + brightness for joining)
"""
import sys
import os
import time
import pickle
import numpy as np
import pandas as pd
import torch
from pathlib import Path

WORK = Path(r"D:\生信\2026Protein Design")
ROUND2 = WORK / "work" / "round2"
ROUND2.mkdir(exist_ok=True, parents=True)

print("="*70)
print("Round 2: ESM2-650M Embeddings for 140K variants")
print("="*70)

# ---- 1. Load GFP_data.xlsx ----
print("\n[1/5] Loading GFP_data.xlsx...", flush=True)
t0 = time.time()
df = pd.read_excel(WORK / "GFP_data.xlsx")
print(f"  Loaded {len(df):,} rows in {time.time()-t0:.1f}s", flush=True)
print(f"  Columns: {list(df.columns)}", flush=True)

# Detect sequence and brightness columns
seq_col = None
bright_col = None
type_col = None
for c in df.columns:
    cl = str(c).lower()
    if 'seq' in cl and seq_col is None:
        seq_col = c
    elif ('bright' in cl or 'fluorescence' in cl or 'signal' in cl) and bright_col is None:
        bright_col = c
    elif ('type' in cl or 'parent' in cl or 'protein' in cl) and type_col is None:
        type_col = c

print(f"  seq_col: {seq_col}", flush=True)
print(f"  bright_col: {bright_col}", flush=True)
print(f"  type_col: {type_col}", flush=True)

# Drop NA sequences
df = df.dropna(subset=[seq_col, bright_col]).reset_index(drop=True)
print(f"  After dropna: {len(df):,} rows", flush=True)

# Save metadata
metadata = df[[seq_col, bright_col]].copy()
if type_col:
    metadata['type'] = df[type_col]
metadata.columns = ['seq', 'brightness'] + (['type'] if type_col else [])
metadata.to_parquet(ROUND2 / "variant_metadata.parquet", index=False)
print(f"  Saved metadata to variant_metadata.parquet", flush=True)

# ---- 2. Load ESM2-650M from local file ----
print("\n[2/5] Loading ESM2-650M (bypassing torch.hub URL check)...", flush=True)
t0 = time.time()
import esm
from esm import ESM2, Alphabet

# Build ESM2 architecture
def build_esm2_650m():
    """Construct ESM2-650M model and load weights from local file."""
    # ESM2-650M architecture (from fair-esm source)
    # num_layers=33, embed_dim=1280, num_heads=20, alphabet_size=33
    alphabet = Alphabet.from_architecture("ESM-1b")
    model = ESM2(
        num_layers=33,
        embed_dim=1280,
        attention_heads=20,
        alphabet=alphabet,
    )
    model_path = r"C:\Users\A\.cache\torch\hub\checkpoints\esm2_t33_650M_UR50D.pt.75c71e41769c4391ba0186bc8c92d0f7.partial"
    sd = torch.load(model_path, map_location="cpu", weights_only=False)
    if isinstance(sd, dict) and 'model' in sd:
        sd = sd['model']
    # Convert from data parallel if needed
    new_sd = {k.replace("module.", ""): v for k, v in sd.items()}
    model.load_state_dict(new_sd, strict=False)
    return model, alphabet

try:
    model, alphabet = build_esm2_650m()
    print(f"  Built model in {time.time()-t0:.1f}s", flush=True)
except Exception as e:
    print(f"  Build failed: {e}", flush=True)
    print(f"  Falling back to esm.pretrained (might be slow)", flush=True)
    t0 = time.time()
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    print(f"  Loaded in {time.time()-t0:.1f}s", flush=True)

print(f"  embed_dim: {model.embed_dim}", flush=True)
print(f"  num_layers: {model.num_layers}", flush=True)
print(f"  params: {sum(p.numel() for p in model.parameters())/1e6:.1f}M", flush=True)

# ---- 3. Move to GPU and test ----
print("\n[3/5] Moving model to GPU...", flush=True)
t0 = time.time()
model = model.eval().cuda()
print(f"  Moved in {time.time()-t0:.1f}s", flush=True)
print(f"  GPU mem: {torch.cuda.memory_allocated()/1e9:.2f} GB", flush=True)

# Quick test
batch_converter = alphabet.get_batch_converter()
test_seq = ('test', 'MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK')
_, _, toks = batch_converter([test_seq])
toks = toks.cuda()
t0 = time.time()
with torch.no_grad():
    out = model(toks, repr_layers=[33], return_contacts=False)
emb = out['representations'][33]
print(f"  Test inference: {time.time()-t0:.3f}s", flush=True)
print(f"  Embedding shape: {emb.shape}", flush=True)

# ---- 4. Batch embed all variants ----
print("\n[4/5] Computing embeddings for all variants...", flush=True)
seqs = df[seq_col].tolist()
print(f"  Total: {len(seqs):,}", flush=True)
batch_size = 32
embed_dim = 1280
all_embs = np.zeros((len(seqs), embed_dim), dtype=np.float32)

t_start = time.time()
for i in range(0, len(seqs), batch_size):
    batch_seqs = seqs[i:i+batch_size]
    batch_data = [(f"p{j}", s) for j, s in enumerate(batch_seqs)]
    try:
        _, _, toks = batch_converter(batch_data)
    except Exception as e:
        # Truncate long sequences
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

    if (i // batch_size) % 50 == 0:
        elapsed = time.time() - t_start
        rate = (i + len(batch_seqs)) / elapsed
        eta = (len(seqs) - i - len(batch_seqs)) / rate
        print(f"  [{i+len(batch_seqs):,}/{len(seqs):,}] {elapsed:.0f}s elapsed, ETA {eta:.0f}s", flush=True)

print(f"\n  Total time: {time.time()-t_start:.0f}s", flush=True)

# ---- 5. Save embeddings ----
print("\n[5/5] Saving embeddings...", flush=True)
np.save(ROUND2 / "esm650m_embeddings.npy", all_embs)
print(f"  Saved {all_embs.shape} to esm650m_embeddings.npy", flush=True)
print(f"  File size: {os.path.getsize(ROUND2 / 'esm650m_embeddings.npy')/1e9:.2f} GB", flush=True)

print("\nDONE", flush=True)