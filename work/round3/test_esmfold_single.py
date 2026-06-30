"""快速测试 ESMFold 单序列折叠"""
import torch
import numpy as np
from transformers import AutoTokenizer, EsmForProteinFolding

seq = "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPKHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

print(f"Loading model...")
tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True)
model = model.cuda()
model.trunk.set_chunk_size(128)
print("Model ready.")

print(f"Folding {len(seq)} aa sequence...")
with torch.no_grad():
    tokens = tokenizer([seq], return_tensors="pt", add_special_tokens=False)["input_ids"].cuda()
    print(f"Tokens shape: {tokens.shape}")
    output = model(tokens)

print(f"Output keys: {list(output.keys())}")

plddt_raw = output["plddt"].cpu().numpy()[0]
atom_mask = output["atom37_atom_exists"].cpu().numpy()[0]
print(f"plddt shape: {plddt_raw.shape}, min: {plddt_raw.min():.4f}, max: {plddt_raw.max():.4f}")
print(f"atom_mask shape: {atom_mask.shape}, num existing: {atom_mask.sum():.0f}/{atom_mask.size}")

# Apply mask and scale
plddt_scaled = plddt_raw * 100.0
masked_sums = (plddt_scaled * atom_mask).sum(axis=1)
masked_counts = atom_mask.sum(axis=1).astype(float)
masked_counts[masked_counts == 0] = 1
plddt_per_res = masked_sums / masked_counts
mean_plddt = float(plddt_per_res.mean())
print(f"Per-residue pLDDT: mean={mean_plddt:.1f}, min={plddt_per_res.min():.1f}, max={plddt_per_res.max():.1f}")

n_low = int((plddt_per_res < 50).sum().item())
n_mid = int(((plddt_per_res >= 50) & (plddt_per_res < 70)).sum().item())
n_high = int(((plddt_per_res >= 70) & (plddt_per_res < 90)).sum().item())
n_vhigh = int((plddt_per_res >= 90).sum().item())
print(f"Distribution: <50:{n_low}  50-70:{n_mid}  70-90:{n_high}  >90:{n_vhigh}")

print("Done!")
