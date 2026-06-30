"""ESMFold pLDDT 验证 (GPU版, FP32, 0-100标度)"""
import json, torch, time
from pathlib import Path
from transformers import AutoTokenizer, EsmForProteinFolding

ROOT = Path(r"D:\生信\2026Protein Design")

with open(ROOT / "work/round3/candidates_round3.json") as f:
    candidates = json.load(f)

print(f"加载 ESMFold 模型 (GPU, FP32)...")
tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained(
    "facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True
).cuda()
# Keep FP32 (FP16 degrades pLDDT significantly)
model.trunk.set_chunk_size(128)
print("模型就绪\n")

results = []
for i, c in enumerate(candidates):
    name = c["name"]
    seq = c["seq"]
    print(f"[{i+1}/{len(candidates)}] Folding {name} ({len(seq)} aa)...", end=" ", flush=True)
    t0 = time.time()

    with torch.no_grad():
        tokens = tokenizer([seq], return_tensors="pt", add_special_tokens=False)["input_ids"].cuda()
        output = model(tokens)

    elapsed = time.time() - t0

    # pLDDT 提取 (shape: [1, L, 37]) 值域 [0,1] -> *100 标度
    plddt_raw = output.plddt.cpu().numpy()[0]
    atom_mask = output.atom37_atom_exists.cpu().numpy()[0]

    plddt_scaled = plddt_raw * 100.0
    masked_sums = (plddt_scaled * atom_mask).sum(axis=1)
    masked_counts = atom_mask.sum(axis=1).astype(float)
    masked_counts[masked_counts == 0] = 1
    plddt_per_res = masked_sums / masked_counts
    mean_plddt = float(plddt_per_res.mean())

    n_low = int((plddt_per_res < 50).sum())
    n_mid = int(((plddt_per_res >= 50) & (plddt_per_res < 70)).sum())
    n_high = int(((plddt_per_res >= 70) & (plddt_per_res < 90)).sum())
    n_vhigh = int((plddt_per_res >= 90).sum())

    status = "✓" if mean_plddt >= 70 else "✗"
    print(f"pLDDT={mean_plddt:.1f} {status}  <50:{n_low}  50-70:{n_mid}  70-90:{n_high}  >90:{n_vhigh}  ({elapsed:.1f}s)")

    results.append({
        "name": name, "length": len(seq),
        "plddt_mean": round(mean_plddt, 1),
        "plddt_lt50": n_low, "plddt_50_70": n_mid,
        "plddt_70_90": n_high, "plddt_gt90": n_vhigh,
        "fold_time_s": round(elapsed, 1)
    })

    del tokens, output
    torch.cuda.empty_cache()

with open(ROOT / "work/round3/esmfold_results.json", "w") as f:
    json.dump(results, f, indent=2, ensure_ascii=False)

print(f"\n结果已保存到 work/round3/esmfold_results.json")
print("\n" + "=" * 60)
print("pLDDT 排序 (mean)")
print("=" * 60)
results.sort(key=lambda x: x["plddt_mean"], reverse=True)
for r in results:
    bar = "█" * int(r["plddt_mean"] / 5)
    print(f"  {r['plddt_mean']:5.1f} {bar} {r['name']}")