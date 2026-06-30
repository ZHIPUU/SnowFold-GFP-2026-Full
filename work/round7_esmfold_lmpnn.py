"""
Round 7: ESMFold 验证 LigandMPNN 候选
"""
import json, time, torch, numpy as np, os
os.environ["CURL_CA_BUNDLE"] = ""
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
R7 = ROOT / "work" / "round7"

with open(R7 / "candidates_lmpnn_r7.json", encoding="utf-8") as f:
    candidates = json.load(f)
print(f"候选: {len(candidates)} 条")

from transformers import AutoTokenizer, EsmForProteinFolding
tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained(
    "facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True
).cuda().eval()

t0 = time.time()
for i, c in enumerate(candidates):
    with torch.no_grad():
        tok = tokenizer([c["seq"]], return_tensors="pt", add_special_tokens=False)
        tok = {k: v.cuda() for k, v in tok.items()}
        out = model(**tok)
    
    plddt = out["plddt"][0, 1:len(c["seq"])+1].cpu().numpy()
    plddt_per_res = plddt.mean(axis=1) * 100
    ptm = out["ptm"].item()
    
    r1, r2 = slice(57, 72), slice(209, 230)
    chromo = np.concatenate([plddt_per_res[r1], plddt_per_res[r2]]).mean()
    
    c["plddt_mean"] = float(plddt_per_res.mean())
    c["plddt_chromo_region"] = float(chromo)
    c["ptm"] = float(ptm)
    c["pass_ptm"] = bool(ptm > 0.75)
    c["pass_plddt"] = bool(plddt_per_res.mean() > 80.0)
    c["pass_chromo"] = bool(chromo > 85.0)
    c["all_pass"] = bool(c["pass_ptm"] and c["pass_plddt"] and c["pass_chromo"])
    
    if (i+1) % 20 == 0:
        print(f"  [{i+1}/{len(candidates)}] {time.time()-t0:.0f}s")

print(f"完成! {time.time()-t0:.0f}s")

candidates.sort(key=lambda x: -(x.get("ptm",0) or 0))
pass_all = [r for r in candidates if r["all_pass"]]

print(f"\n{'='*80}")
print(f"结果: 全部通过={len(pass_all)}, 仅pTM={sum(1 for r in candidates if r['pass_ptm'])}/{len(candidates)}")
print(f"{'='*80}")

print(f"\n{'#':<4} {'Name':<28} {'pLDDT':>6} {'Chromo':>7} {'pTM':>6} {'Mut':>4} {'Status'}")
for i, r in enumerate(candidates[:20], 1):
    ptm = r.get("ptm",0) or 0
    plddt = r.get("plddt_mean",0) or 0
    chromo = r.get("plddt_chromo_region",0) or 0
    if r["all_pass"]: s = "✅ ALL"
    elif r["pass_ptm"]: s = "🟡 pTM"
    elif ptm > 0.7: s = "🟠 Close"
    else: s = "🔴 No"
    print(f"{i:<4} {r['name'][:28]:<28} {plddt:>6.1f} {chromo:>7.1f} {ptm:>6.4f} {r['n_muts']:>4} {s}")

with open(R7 / "esmfold_lmpnn_r7.json", "w", encoding="utf-8") as f:
    json.dump(candidates, f, indent=2, ensure_ascii=False)

# 对比 ProteinMPNN 最佳
print(f"\n对比 MPNN_T01_014 (ProteinMPNN 最佳):")
print(f"  pLDDT=68.3, chromo=67.6, pTM=0.765")
if candidates:
    print(f"  LMPNN 最佳: pLDDT={candidates[0]['plddt_mean']:.1f}, chromo={candidates[0]['plddt_chromo_region']:.1f}, pTM={candidates[0]['ptm']:.4f}")
