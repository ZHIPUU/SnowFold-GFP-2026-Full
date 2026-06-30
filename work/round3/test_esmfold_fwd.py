"""Minimal ESMFold forward test"""
import torch
import time
from transformers import AutoTokenizer, EsmForProteinFolding

seq = "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPKHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

print("Loading model...")
t0 = time.time()
tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True)
model = model.cuda()
# DO NOT set chunk_size - use default
t1 = time.time()
print(f"Model loaded in {t1-t0:.1f}s")

print(f"Folding {len(seq)} aa sequence...")
with torch.no_grad():
    tokens = tokenizer([seq], return_tensors="pt", add_special_tokens=False)["input_ids"].cuda()
    print(f"Tokens shape: {tokens.shape}")
    t2 = time.time()
    output = model(tokens)
    t3 = time.time()
    print(f"Forward pass took {t3-t2:.1f}s")

print(f"Output keys: {list(output.keys())}")
plddt = output["plddt"]
print(f"pLDDT shape: {plddt.shape}, values: min={plddt.min():.4f}, max={plddt.max():.4f}, mean={plddt.mean():.4f}")
print("Done!")
