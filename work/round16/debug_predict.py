"""Debug predict_and_save"""
import warnings, time, torch, torch.nn.functional as F, numpy as np
warnings.filterwarnings("ignore")
from transformers import AutoTokenizer, EsmForProteinFolding
import json

# Load
tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True).cuda()
model.trunk.set_chunk_size(128); model.eval()

# Load seq
data = json.load(open(r"D:\生信\2026Protein Design\work\round15\final_6_r15.json"))
seq = data[0]["seq"]
print(f"Seq len={len(seq)}, starts with M: {seq.startswith('M')}")

# ESMFold
inputs = tok([seq], return_tensors="pt", add_special_tokens=False)
inputs = {k: v.cuda() for k, v in inputs.items()}
print("Running ESMFold...", flush=True)
t0 = time.time()
with torch.no_grad():
    out = model(**inputs, num_recycles=8)
print(f"ESMFold: {time.time()-t0:.1f}s", flush=True)

# pLDDT
logits = out.plddt[0]
probs = F.softmax(logits, dim=-1)
centers = torch.linspace(0.5/37, 1-0.5/37, 37, device=logits.device)
plddt_01 = (probs * centers.unsqueeze(0)).sum(-1)
plddt_100 = plddt_01.cpu().numpy() * 100
print(f"pLDDT: global={plddt_100.mean():.1f}, chromo[57:72]={plddt_100[57:72].mean():.1f}, pTM={out.ptm.item():.4f}")

# Positions
positions = out.positions[-1][0].cpu().numpy()
print(f"Positions shape: {positions.shape}")

# Save PDB
aa3 = {a:b for a,b in zip("ACDEFGHIKLMNPQRSTVWY", "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split())}
pdb_path = r"D:\生信\2026Protein Design\work\round16\pdbs\debug_test.pdb"
with open(pdb_path, "w") as f:
    f.write("REMARK  debug\n")
    aidx = 1
    for i in range(len(seq)):
        rn = aa3.get(seq[i], 'ALA')
        for j, an in enumerate(["N","CA","C","O"]):
            x, y, z = positions[i, j]
            f.write(f"ATOM  {aidx:5d} {an:^4s} {rn:3s} A{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{plddt_100[i]:6.2f}\n")
            aidx += 1
    f.write("END\n")
print(f"PDB saved: {pdb_path}")
print("ALL OK")
