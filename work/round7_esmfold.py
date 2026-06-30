"""Round 7: ESMFold 批量验证 (HuggingFace API)"""
import json, time, torch, numpy as np
from pathlib import Path
import os
os.environ["CURL_CA_BUNDLE"] = ""

ROOT = Path(r"D:\生信\2026Protein Design")
R7 = ROOT / "work" / "round7"

with open(R7 / "candidates_round7.json", encoding="utf-8") as f:
    candidates = json.load(f)
print(f"候选: {len(candidates)} 条")

print("\n加载 ESMFold (HuggingFace)...")
from transformers import AutoTokenizer, EsmForProteinFolding

tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained(
    "facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True
).cuda().eval()
print("OK")

batch_size = 10
t0 = time.time()

for i, c in enumerate(candidates):
    seq = c["seq"]
    
    with torch.no_grad():
        tokenized = tokenizer([seq], return_tensors="pt", add_special_tokens=False)
        tokenized = {k: v.cuda() for k, v in tokenized.items()}
        output = model(**tokenized)
    
    # pLDDT: HF EsmForProteinFolding 返回的 plddt 是 0-1 范围的期望值
    # 需要乘以 100 得到标准 0-100 pLDDT
    plddt_raw = output["plddt"][0, 1:len(seq)+1].cpu().numpy()  # (L, 37) per-atom, 0-1
    plddt_per_res = plddt_raw.mean(axis=1) * 100  # 0-100 scale
    ptm = output["ptm"].item()  # pTM 已经是 0-1 范围
    
    # 生色团区域
    region1 = slice(57, 72)
    region2 = slice(209, 230)
    chromo = np.concatenate([plddt_per_res[region1], plddt_per_res[region2]]).mean() if len(plddt_per_res) > 230 else float(plddt_per_res.mean())
    
    c["plddt_mean"] = float(plddt_per_res.mean())
    c["plddt_chromo_region"] = float(chromo)
    c["ptm"] = float(ptm)
    c["pass_ptm"] = bool(ptm > 0.75)
    c["pass_plddt"] = bool(plddt_per_res.mean() > 80.0)
    c["pass_chromo"] = bool(chromo > 85.0)
    c["all_pass"] = c["pass_ptm"] and c["pass_plddt"] and c["pass_chromo"]
    
    if (i+1) % batch_size == 0:
        elapsed = time.time() - t0
        print(f"  [{i+1}/{len(candidates)}] {elapsed:.0f}s ({elapsed/(i+1):.1f}s/seq)")

print(f"\n完成! {time.time()-t0:.0f}s")

candidates.sort(key=lambda x: -(x.get("ptm",0) or 0))
pass_all = [r for r in candidates if r["all_pass"]]

print(f"\n{'='*80}")
print(f"结果: 全部通过={len(pass_all)}, 仅pTM={sum(1 for r in candidates if r['pass_ptm'])} / {len(candidates)}")
print(f"{'='*80}")

print(f"\n{'#':<4} {'Name':<25} {'pLDDT':>6} {'Chromo':>7} {'pTM':>6} {'Mut':>4} {'Status'}")
for i, r in enumerate(candidates[:20], 1):
    ptm = r.get("ptm",0) or 0
    plddt = r.get("plddt_mean",0) or 0
    chromo = r.get("plddt_chromo_region",0) or 0
    
    if r["all_pass"]: s = "✅ ALL"
    elif r["pass_ptm"]: s = "🟡 pTM"
    elif ptm > 0.7: s = "🟠 Close"
    else: s = "🔴 No"
    
    print(f"{i:<4} {r['name'][:25]:<25} {plddt:>6.1f} {chromo:>7.1f} {ptm:>6.4f} {r['n_muts']:>4} {s}")

with open(R7 / "esmfold_round7.json", "w", encoding="utf-8") as f:
    json.dump(candidates, f, indent=2, ensure_ascii=False)
print(f"\n已保存: {R7 / 'esmfold_round7.json'}")
