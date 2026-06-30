"""Quick ESMFold test"""
import warnings, torch, time
warnings.filterwarnings("ignore")
from transformers import AutoTokenizer, EsmForProteinFolding

print("Loading ESMFold...")
tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True).cuda()
model.trunk.set_chunk_size(128)
model.eval()
print("Loaded.")

seq = "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
inputs = tok([seq], return_tensors="pt", add_special_tokens=False)
inputs = {k: v.cuda() for k, v in inputs.items()}
print(f"Input: seq len={len(seq)}, input_ids shape={inputs['input_ids'].shape}")

t0 = time.time()
with torch.no_grad():
    out = model(**inputs, num_recycles=8)
t1 = time.time()
print(f"ESMFold took {t1-t0:.1f}s")
print(f"pLDDT shape: {out.plddt.shape}")
print(f"pTM: {out.ptm.item():.4f}")
print(f"positions: {len(out.positions[-1][0])} residues")
print(f"positions[0] shape: {out.positions[-1][0][0].shape}")

plddt = out.plddt[0, :, 1].cpu().numpy()
print(f"Global pLDDT mean: {plddt.mean():.1f}")
print(f"Chromo region (57:72) pLDDT: {plddt[57:72].mean():.1f}")
print(f"SUCCESS - ESMFold works fine")
