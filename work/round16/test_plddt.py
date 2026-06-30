"""Verify correct pLDDT computation from ESMFold"""
import warnings, torch, torch.nn.functional as F
warnings.filterwarnings("ignore")
from transformers import AutoTokenizer, EsmForProteinFolding

tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True).cuda()
model.trunk.set_chunk_size(128)
model.eval()

# sfGFP WT - expected pLDDT ~48-50
seq = "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
inputs = tok([seq], return_tensors="pt", add_special_tokens=False)
inputs = {k: v.cuda() for k, v in inputs.items()}
with torch.no_grad():
    out = model(**inputs, num_recycles=8)

# Correct pLDDT computation
logits = out.plddt[0]  # (L, 37)
probs = F.softmax(logits, dim=-1)
bin_width = 1.0 / 37
bin_centers = torch.linspace(0.5 * bin_width, 1.0 - 0.5 * bin_width, 37, device=logits.device)
plddt_01 = (probs * bin_centers.unsqueeze(0)).sum(dim=-1)  # (L,) 0-1
plddt_100 = plddt_01.cpu().numpy() * 100

ptm = float(out.ptm.cpu().item())
print(f"sfGFP WT: Global pLDDT={plddt_100.mean():.1f}, Chromo[57:72]={plddt_100[57:72].mean():.1f}, pTM={ptm:.3f}")
print(f"Expected: Global pLDDT ~48-50")
print(f"PASS: Score computation correct" if 40 < plddt_100.mean() < 60 else "WARNING: Unexpected values")

# Also test with the R15 Top 1 sequence
r15_seq = "MTIPGETLLAGVVPVRVNLDGDVNGKKFKVVGEGEGDATKGELRLTFEVTEGELPLDPLLLSYILTYPLSIFRKMPKDNPLRPFYLACLPEGYVIERTLDFKGEGTLNVTSETYFDGDTLVSNITLKGTGFKEGGKLLTKKVKEIRLVGDLTITPDEEKKGVKLTYTLELTFEDGSTSTADVKELIYPKGKGPETLPEPQTLTVDLTLTAVPEVEGDRFRFEHRDPPGVEPPPLSELSK"
inputs2 = tok([r15_seq], return_tensors="pt", add_special_tokens=False)
inputs2 = {k: v.cuda() for k, v in inputs2.items()}
with torch.no_grad():
    out2 = model(**inputs2, num_recycles=8)

logits2 = out2.plddt[0]
probs2 = F.softmax(logits2, dim=-1)
plddt2_01 = (probs2 * bin_centers.unsqueeze(0)).sum(dim=-1)
plddt2_100 = plddt2_01.cpu().numpy() * 100
ptm2 = float(out2.ptm.cpu().item())
chromo2 = plddt2_100[57:72].mean()
score = 0.40 * ptm2 + 0.30 * (plddt2_100.mean() / 100) + 0.30 * (chromo2 / 100)
print(f"\nR15 Top 1: Global pLDDT={plddt2_100.mean():.1f}, Chromo={chromo2:.1f}, pTM={ptm2:.3f}, Score={score:.4f}")
print(f"Expected: pLDDT ~84, Chromo ~89, pTM ~0.90, Score ~0.88")
