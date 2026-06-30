import torch
import esm
import time

print("Loading ESM2-150M...")
t0 = time.time()
model, alphabet = esm.pretrained.esm2_t30_150M_UR50D()
batch_converter = alphabet.get_batch_converter()
print(f"Load time: {time.time()-t0:.1f}s")

device = "cuda"
model = model.to(device).eval()
print(f"Model on {device}")

seq = ("test", "MASKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK")
batch_labels, batch_strs, batch_tokens = batch_converter([seq])
batch_tokens = batch_tokens.to(device)
torch.cuda.synchronize()
t0 = time.time()
with torch.no_grad():
    results = model(batch_tokens, repr_layers=[30], return_contacts=False)
torch.cuda.synchronize()
print(f"Inference time: {time.time()-t0:.2f}s")
emb = results["representations"][30]
print(f"Embeddings shape: {emb.shape}")
print(f"GPU mem: {torch.cuda.memory_allocated()/1024**3:.2f} GB")
print("ESM2-150M OK!")

# Test batching 4 sequences
print("\nTesting batch=4...")
seqs = [("s1", seq[1])] * 4
batch_labels, batch_strs, batch_tokens = batch_converter(seqs)
batch_tokens = batch_tokens.to(device)
torch.cuda.synchronize()
t0 = time.time()
with torch.no_grad():
    results = model(batch_tokens, repr_layers=[30], return_contacts=False)
torch.cuda.synchronize()
print(f"Batch=4 inference: {time.time()-t0:.2f}s")
print(f"GPU mem: {torch.cuda.memory_allocated()/1024**3:.2f} GB")
