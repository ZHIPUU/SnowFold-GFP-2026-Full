"""Verify ESM2-650M model file integrity."""
import os
import hashlib
import time

p = r'C:\Users\A\.cache\torch\hub\checkpoints\esm2_t33_650M_UR50D.pt.75c71e41769c4391ba0186bc8c92d0f7.partial'

print(f"File: {p}")
print(f"Size: {os.path.getsize(p)/1e9:.3f} GB ({os.path.getsize(p):,} bytes)")
print()

print("Computing SHA256 (this takes a while)...")
t0 = time.time()
h = hashlib.sha256()
with open(p, 'rb') as f:
    for chunk in iter(lambda: f.read(1024*1024*32), b''):
        h.update(chunk)
print(f"SHA256: {h.hexdigest()}")
print(f"Time: {time.time()-t0:.1f}s")