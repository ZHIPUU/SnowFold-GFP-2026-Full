"""Test load ESM2-650M with progress tracking."""
import sys
import time

print("Step 1: Importing esm...", flush=True)
t0 = time.time()
import esm
print(f"  Done in {time.time()-t0:.1f}s", flush=True)

print("Step 2: Loading esm2_t33_650M_UR50D model...", flush=True)
t0 = time.time()
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
print(f"  Loaded in {time.time()-t0:.1f}s", flush=True)
print(f"  embed_dim: {model.embed_dim}", flush=True)
print(f"  layers: {model.num_layers}", flush=True)
print(f"  params: {sum(p.numel() for p in model.parameters())/1e6:.1f}M", flush=True)

print("Step 3: Moving to GPU...", flush=True)
t0 = time.time()
model = model.eval().cuda()
print(f"  Moved in {time.time()-t0:.1f}s", flush=True)

print("Step 4: Test inference...", flush=True)
batch_converter = alphabet.get_batch_converter()
seq = ('test', 'MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK')
labels, strs, toks = batch_converter([seq])
toks = toks.cuda()
t0 = time.time()
import torch
with torch.no_grad():
    out = model(toks, repr_layers=[33], return_contacts=False)
print(f"  Inference: {time.time()-t0:.3f}s", flush=True)
emb = out['representations'][33]
print(f"  Embedding shape: {emb.shape}", flush=True)
print(f"  GPU mem: {torch.cuda.memory_allocated()/1e9:.2f} GB", flush=True)
print("ALL OK", flush=True)