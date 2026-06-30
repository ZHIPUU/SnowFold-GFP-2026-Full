"""Debug ESMFold pLDDT output structure"""
import torch, numpy as np
from transformers import AutoTokenizer, EsmForProteinFolding
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")

with open(ROOT / "AAseqs of 5 GFP proteins_20260511.txt") as f:
    wt_text = f.read()

avGFP = ""
for block in wt_text.split(">"):
    block = block.strip()
    if not block or "avGFP" not in block:
        continue
    lines = block.split("\n")
    avGFP = "".join(l.strip() for l in lines[1:] if l.strip() and not l.startswith("#"))

print(f"avGFP: {len(avGFP)} aa")

print("Loading ESMFold...")
tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
model = EsmForProteinFolding.from_pretrained(
    "facebook/esmfold_v1", low_cpu_mem_usage=True
).cuda()
model.esm = model.esm.half()
model.trunk.set_chunk_size(64)
print("Ready\n")

with torch.no_grad():
    tokens = tokenizer([avGFP], return_tensors="pt", add_special_tokens=False)["input_ids"].cuda()
    output = model(tokens)

# Debug
print(f"Output type: {type(output).__name__}")
print(f"Output attributes: {[k for k in dir(output) if not k.startswith('_')]}\n")

# Check plddt via attribute
if hasattr(output, "plddt"):
    x = output.plddt
    print(f"plddt type: {type(x)}, shape: {x.shape}")
    xnp = x.cpu().numpy()
    print(f"  numpy shape: {xnp.shape}")
    print(f"  mean: {xnp.mean():.3f}")
    # Take first dimension if it exists
    if xnp.ndim >= 2:
        xnp = xnp[0]
    print(f"  per-residue (first 10): {xnp[:10]}")
    print(f"  per-residue mean: {xnp.mean():.3f}")
    print(f"  min: {xnp.min():.3f}, max: {xnp.max():.3f}")

# Check via dict
try:
    x2 = output["plddt"]
    print(f"\nplddt via dict: type={type(x2)}, shape={x2.shape}")
except Exception as e:
    print(f"\nplddt via dict: ERROR - {e}")

# Check for atom37
if hasattr(output, "positions"):
    print(f"\npositions shape: {output.positions.shape}")

# Check ptm
if hasattr(output, "ptm"):
    ptm = output.ptm
    print(f"\nptm: {ptm}")

# Try to find the correct pLDDT
if hasattr(output, "atom37_lddt"):
    print(f"\natom37_lddt shape: {output.atom37_lddt.shape}")

# Let's look at any attribute containing 'lddt' or 'plddt'
for attr in dir(output):
    if not attr.startswith("_") and ("lddt" in attr.lower() or "plddt" in attr.lower()):
        val = getattr(output, attr)
        if hasattr(val, "shape"):
            print(f"\n  {attr}: shape={val.shape}, dtype={val.dtype}")
