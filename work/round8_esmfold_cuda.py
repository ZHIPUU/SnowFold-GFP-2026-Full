"""
Round 8: ESMFold验证 LigandMPNN 候选 (CUDA)
"""
import json, time, torch, numpy as np
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
R8 = ROOT / "work" / "round8"

with open(R8 / "candidates_lmpnn_r8.json", encoding="utf-8") as f:
    candidates = json.load(f)
print(f"Candidates: {len(candidates)}")

from transformers import AutoTokenizer, EsmForProteinFolding

tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained(
    "facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True
).cuda().eval()
print(f"ESMFold on CUDA: {next(model.parameters()).is_cuda}")

t0 = time.time()
for i, c in enumerate(candidates):
    with torch.no_grad():
        tok = tokenizer([c["seq"]], return_tensors="pt", add_special_tokens=False)
        tok = {k: v.cuda() for k, v in tok.items()}
        out = model(**tok)
    
    plddt = out["plddt"][0, 1:len(c["seq"])+1].cpu().numpy()
    plddt_r = plddt.mean(axis=1) * 100
    ptm = out["ptm"].item()
    r1, r2 = slice(57, 72), slice(209, 230)
    chromo = np.concatenate([plddt_r[r1], plddt_r[r2]]).mean() if len(plddt_r) > 230 else float(plddt_r.mean())
    
    c["plddt_mean"] = float(plddt_r.mean())
    c["plddt_chromo_region"] = float(chromo)
    c["ptm"] = float(ptm)
    c["pass_ptm"] = bool(ptm > 0.75)
    c["pass_plddt"] = bool(plddt_r.mean() > 80.0)
    c["pass_chromo"] = bool(chromo > 85.0)
    
    if (i+1) % 30 == 0:
        print(f"  [{i+1}/{len(candidates)}] {time.time()-t0:.0f}s")

print(f"Done! {time.time()-t0:.0f}s")

# Sort by sort_score (new rules)
for c in candidates:
    ptm = c.get("ptm", 0) or 0
    plddt = c.get("plddt_mean", 0) or 0
    chromo = c.get("plddt_chromo_region", plddt) or plddt
    nmuts = c.get("n_muts", 0) or 0
    c["sort_score"] = round(ptm * 0.50 + plddt / 100 * 0.30 + chromo / 100 * 0.20 - nmuts / 200 * 0.05, 4)

candidates.sort(key=lambda x: -x["sort_score"])

top2 = [c for c in candidates if (c.get("ptm") or 0) > 0.70 and (c.get("plddt_mean") or 0) > 60]
remaining = [c for c in candidates if (c.get("ptm") or 0) > 0.55]

print(f"\nTop2达标: {len(top2)}")
print(f"剩余达标: {len(remaining)}")
print(f"\n{'#':<4} {'Name':<30} {'pLDDT':>6} {'Chromo':>6} {'pTM':>6} {'Muts':>4} {'Score':>8}")
for i, c in enumerate(candidates[:20], 1):
    print(f"{i:<4} {c['name'][:30]:<30} {c['plddt_mean']:>6.1f} {c.get('plddt_chromo_region',0):>6.1f} {c['ptm']:>6.4f} {c['n_muts']:>4} {c['sort_score']:>8.4f}")

# Compare with Round 4 best
print(f"\n对比 R4: MPNN_T01_014 (pLDDT=68.3, chromo=67.6, pTM=0.765)")
if candidates:
    print(f"R8 LMPNN最佳: pLDDT={candidates[0]['plddt_mean']:.1f}, chromo={candidates[0].get('plddt_chromo_region',0):.1f}, pTM={candidates[0]['ptm']:.4f}")

with open(R8 / "esmfold_lmpnn_r8.json", "w", encoding="utf-8") as f:
    json.dump(candidates, f, indent=2, ensure_ascii=False)
print(f"\nSaved: {R8 / 'esmfold_lmpnn_r8.json'}")
