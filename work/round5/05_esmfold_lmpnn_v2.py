"""Round 5 P0-1d: LigandMPNN v2 候选 ESMFold 评估"""
import json, torch, time
from pathlib import Path
from transformers import AutoTokenizer, EsmForProteinFolding

ROOT = Path(r"D:\生信\2026Protein Design")
R5 = ROOT / "work" / "round5"

with open(R5 / "lmpnn_v2_candidates.json", encoding="utf-8") as f:
    cands = json.load(f)
print(f"候选: {len(cands)}")

print(f"\nGPU: {torch.cuda.get_device_name(0)}")
tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained(
    "facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True).cuda()
model.trunk.set_chunk_size(128); model.eval()
print("模型就绪\n")

results = []
for i, c in enumerate(cands):
    seq = c["seq"]
    print(f"[{i+1}/{len(cands)}] {c['name']} ({len(seq)} aa)...", end=" ", flush=True)
    t0 = time.time()
    with torch.no_grad():
        tokens = tokenizer([seq], return_tensors="pt", add_special_tokens=False)["input_ids"].cuda()
        output = model(tokens)
    plddt_raw = output.plddt.cpu().numpy()[0]
    atom_mask = output.atom37_atom_exists.cpu().numpy()[0]
    plddt_scaled = plddt_raw * 100.0
    masked_sums = (plddt_scaled * atom_mask).sum(axis=1)
    masked_counts = atom_mask.sum(axis=1).astype(float)
    masked_counts[masked_counts == 0] = 1
    plddt_per_res = masked_sums / masked_counts
    mean_plddt = float(plddt_per_res.mean())
    cb_region = plddt_per_res[58:73]
    cb_mean = float(cb_region.mean())
    try:
        ptm = float(output.ptm.cpu().item())
    except:
        ptm = None
    elapsed = time.time() - t0
    ptm_s = f"{ptm:.3f}" if ptm else "N/A"
    print(f"pLDDT={mean_plddt:5.1f} cb={cb_mean:5.1f} pTM={ptm_s} ({elapsed:.1f}s)")
    c["plddt_mean"] = round(mean_plddt, 2)
    c["plddt_chromo_region"] = round(cb_mean, 2)
    c["ptm"] = round(ptm, 4) if ptm else None
    c["fold_time_s"] = round(elapsed, 1)
    results.append(c)
    del tokens, output
    torch.cuda.empty_cache()

with open(R5 / "esmfold_lmpnn_v2.json", "w", encoding="utf-8") as f:
    json.dump(results, f, indent=2, ensure_ascii=False)

print(f"\n保存 esmfold_lmpnn_v2.json ({len(results)} 条)")
print(f"\nTop 10 by pLDDT:")
for r in sorted(results, key=lambda x: -x["plddt_mean"])[:10]:
    print(f"  pLDDT={r['plddt_mean']:5.1f} cb={r['plddt_chromo_region']:5.1f} pTM={r['ptm']:.3f}  {r['name']} (n_mut={r['n_muts']}, lig={r['lmpnn_ligand']:.3f})")
