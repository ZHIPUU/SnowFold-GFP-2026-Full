"""
Round 4 Step 2B: GPU ESMFold 批量评估扩展候选池
==================================
评估全部 22 条候选 (修正后的)
"""
import json, torch, time
from pathlib import Path
from transformers import AutoTokenizer, EsmForProteinFolding

ROOT = Path(r"D:\生信\2026Protein Design")
OUT = ROOT / "work" / "round4"

with open(OUT / "candidates_round4_extended.json", encoding="utf-8") as f:
    candidates = json.load(f)

print(f"扩展候选池: {len(candidates)} 条")
from collections import Counter
for role, count in Counter(c["role"] for c in candidates).most_common():
    print(f"  {role}: {count}")

# ============================================================
print(f"\nGPU: {torch.cuda.get_device_name(0)}")
print("加载 ESMFold (FP32)...")
tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained(
    "facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True
).cuda()
model.trunk.set_chunk_size(128)
model.eval()
print("✓ 模型就绪\n")

# ============================================================
results = []
total_start = time.time()
for i, c in enumerate(candidates):
    name = c["name"]
    seq = c["seq"]
    print(f"[{i+1:>2}/{len(candidates)}] {name:<32} ({len(seq)} aa)...", end=" ", flush=True)
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

    # chromophore 区域: sfGFP/avGFP 是 65-67, mBaoJin 是 66-68
    if c["scaffold"] == "mBaoJin":
        cb_start, cb_end = 60, 75  # 包含 GYG + 邻位
    else:
        cb_start, cb_end = 58, 73  # TYG + 邻位
    cb_region = plddt_per_res[max(0, cb_start):min(len(plddt_per_res), cb_end)]
    cb_mean = float(cb_region.mean())

    n_low = int((plddt_per_res < 50).sum())
    n_mid = int(((plddt_per_res >= 50) & (plddt_per_res < 70)).sum())
    n_high = int(((plddt_per_res >= 70) & (plddt_per_res < 90)).sum())
    n_vhigh = int((plddt_per_res >= 90).sum())

    try:
        ptm = float(output.ptm.cpu().item()) if hasattr(output, 'ptm') else None
    except Exception:
        ptm = None

    elapsed = time.time() - t0
    icon = "🟢" if mean_plddt >= 70 else ("🟡" if mean_plddt >= 50 else "🔴")
    ptm_str = f"{ptm:.3f}" if ptm is not None else "N/A"
    print(f"pLDDT={mean_plddt:5.1f} cb={cb_mean:5.1f} pTM={ptm_str} {icon} ({elapsed:.1f}s)")

    results.append({
        "name": name,
        "role": c["role"],
        "scaffold": c["scaffold"],
        "n_muts": c["n_muts"],
        "length": len(seq),
        "seq": seq,
        "notes": c.get("notes", ""),
        "expected_tm": c.get("expected_tm", 0),
        "plddt_mean": round(mean_plddt, 2),
        "plddt_chromo_region": round(cb_mean, 2),
        "plddt_lt50": n_low,
        "plddt_50_70": n_mid,
        "plddt_70_90": n_high,
        "plddt_gt90": n_vhigh,
        "ptm": round(ptm, 4) if ptm is not None else None,
        "fold_time_s": round(elapsed, 1),
    })

    del tokens, output
    torch.cuda.empty_cache()

print(f"\n总耗时: {(time.time()-total_start)/60:.1f} min")

with open(OUT / "esmfold_round4_extended.json", "w", encoding="utf-8") as f:
    json.dump(results, f, indent=2, ensure_ascii=False)

# ============================================================
# 排序展示
print("\n" + "=" * 100)
print(f"{'name':<32} {'role':<22} {'mut':>3} {'pLDDT':>6} {'cb':>6} {'pTM':>6}")
print("=" * 100)
for r in sorted(results, key=lambda x: -x["plddt_mean"]):
    icon = "🟢" if r["plddt_mean"] >= 70 else ("🟡" if r["plddt_mean"] >= 50 else "🔴")
    ptm_str = f"{r['ptm']:.3f}" if r['ptm'] is not None else "  N/A"
    print(f"{r['name']:<32} {r['role'][:22]:<22} {r['n_muts']:>3} {r['plddt_mean']:>6.2f} {r['plddt_chromo_region']:>6.2f} {ptm_str:>6} {icon}")

print("\n✓ 结果保存到 esmfold_round4_extended.json")
print("下一步: 重跑 03b_score_and_select 选 Top-6")
