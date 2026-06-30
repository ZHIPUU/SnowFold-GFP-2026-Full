"""Load ESM2-650M model directly from local file (no URL check)."""
import sys
import time

print("Step 1: Importing...", flush=True)
t0 = time.time()
import torch
print(f"  torch: {time.time()-t0:.1f}s", flush=True)

t0 = time.time()
from esm import pretrained, Alphabet, ESM2
print(f"  esm: {time.time()-t0:.1f}s", flush=True)

print("Step 2: Loading state dict from local file...", flush=True)
t0 = time.time()
model_path = r"C:\Users\A\.cache\torch\hub\checkpoints\esm2_t33_650M_UR50D.pt.75c71e41769c4391ba0186bc8c92d0f7.partial"
print(f"  File size: {torch.load.__module__}", flush=True)
# Rename to expected name first
import os, shutil
expected = r"C:\Users\A\.cache\torch\hub\checkpoints\esm2_t33_650M_UR50D.pt"
if not os.path.exists(expected):
    print(f"  Renaming partial to expected name...", flush=True)
    shutil.copy(model_path, expected)
    print(f"  Copied to {expected}", flush=True)

print("Step 3: torch.load state dict...", flush=True)
t0 = time.time()
model_data = torch.load(expected, map_location="cpu", weights_only=False)
print(f"  Loaded state dict in {time.time()-t0:.1f}s", flush=True)
print(f"  Type: {type(model_data)}", flush=True)
if isinstance(model_data, tuple):
    print(f"  Tuple length: {len(model_data)}", flush=True)
    for i, x in enumerate(model_data):
        if isinstance(x, dict):
            print(f"    [{i}] dict with {len(x)} keys", flush=True)
        elif isinstance(x, torch.Tensor):
            print(f"    [{i}] tensor {x.shape}", flush=True)
        else:
            print(f"    [{i}] {type(x)}", flush=True)
elif isinstance(model_data, dict):
    print(f"  Keys: {list(model_data.keys())[:5]}...", flush=True)
    print(f"  Total params: {sum(v.numel() if isinstance(v, torch.Tensor) else 0 for v in model_data.values())/1e6:.1f}M", flush=True)

print("DONE", flush=True)