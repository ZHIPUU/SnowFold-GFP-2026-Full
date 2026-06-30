"""ESM pipeline smoke test - 用最小模型验证"""
import time
import torch
import esm

print("Loading ESM2-8M (tiny, ~30MB)...")
t0 = time.time()
model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
batch_converter = alphabet.get_batch_converter()
print(f"  Load: {time.time()-t0:.1f}s")

device = "cuda"
model = model.to(device).eval()

# 加载4个 GFP WT 序列做测试
seqs = [
    ("avGFP", "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"),
    ("amacGFP", "MSKGEELFTGIVPVLIELDGDVHGHKFSVRGEGEGDADYGKLEIKFICTTGKLPVPWPTLVTTLSYGILCFARYPEHMKMNDFFKSAMPEGYIQERTIFFQDDGKYKTRGEVKFEGDTLVNRIELKGMDFKEDGNILGHKLEYNFNSHNVYIMPDKANNGLKVNFKIRHNIEGGGVQLADHYQTNVPLGDGPVLIPINHYLSCQTAISKDRNETRDHMVFLEFFSACGHTHGMDELYK"),
]
data = [(s[0], s[1][:230]) for s in seqs]  # 截断到 230aa
batch_labels, batch_strs, batch_tokens = batch_converter(data)
batch_tokens = batch_tokens.to(device)

torch.cuda.synchronize()
t0 = time.time()
with torch.no_grad():
    results = model(batch_tokens, repr_layers=[6], return_contacts=False)
torch.cuda.synchronize()
print(f"  Inference: {time.time()-t0:.2f}s")
print(f"  Output shape: {results['representations'][6].shape}")
print(f"  GPU mem: {torch.cuda.memory_allocated()/1024**3:.2f} GB")
print("ESM2-8M pipeline OK!")
print(f"\n=== Going for ESM2-150M ===")
del model
torch.cuda.empty_cache()

t0 = time.time()
model, alphabet = esm.pretrained.esm2_t30_150M_UR50D()
batch_converter = alphabet.get_batch_converter()
print(f"  Load: {time.time()-t0:.1f}s")

model = model.to(device).eval()

torch.cuda.synchronize()
t0 = time.time()
with torch.no_grad():
    results = model(batch_tokens, repr_layers=[30], return_contacts=False)
torch.cuda.synchronize()
print(f"  Inference: {time.time()-t0:.2f}s")
print(f"  GPU mem: {torch.cuda.memory_allocated()/1024**3:.2f} GB")
print("ESM2-150M OK!")