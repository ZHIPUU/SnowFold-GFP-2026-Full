"""ESMFold pLDDT 验证 (CPU版)"""
import json
import torch
import numpy as np
from pathlib import Path
from transformers import AutoTokenizer, EsmForProteinFolding
import time

ROOT = Path(r"D:\生信\2026Protein Design")

# 加载候选
with open(ROOT / "work/round3/candidates_round3.json") as f:
    candidates = json.load(f)

print(f"加载 ESMFold 模型 (CPU模式)...")
tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained(
    "facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True
)
# 保持CPU (不要 .cuda())
model.esm = model.esm.float()  # FP32确保精度
model.trunk.set_chunk_size(64)  # CPU下小一点避免内存爆
print("模型就绪\n")

results = []
for i, c in enumerate(candidates):
    name = c["name"]
    seq = c["seq"]
    print(f"[{i+1}/{len(candidates)}] Folding {name} ({len(seq)} aa)...", end=" ", flush=True)
    t0 = time.time()

    with torch.no_grad():
        tokens = tokenizer([seq], return_tensors="pt", add_special_tokens=False)["input_ids"]
        output = model(tokens)

    elapsed = time.time() - t0

    # 提取 pLDDT (shape: [1, L, 37]) 值域 [0,1]
    plddt_raw = output["plddt"].numpy()[0]
    atom_mask = output["atom37_atom_exists"].numpy()[0]

    plddt_scaled = plddt_raw * 100.0
    masked_sums = (plddt_scaled * atom_mask).sum(axis=1)
    masked_counts = atom_mask.sum(axis=1).astype(float)
    masked_counts[masked_counts == 0] = 1
    plddt_per_res = masked_sums / masked_counts
    mean_plddt = float(plddt_per_res.mean())

    n_low = int((plddt_per_res < 50).sum().item())
    n_mid = int(((plddt_per_res >= 50) & (plddt_per_res < 70)).sum().item())
    n_high = int(((plddt_per_res >= 70) & (plddt_per_res < 90)).sum().item())
    n_vhigh = int((plddt_per_res >= 90).sum().item())

    status = "✓" if mean_plddt >= 70 else "✗"
    print(f"pLDDT={mean_plddt:.1f} {status}  <50:{n_low}  50-70:{n_mid}  70-90:{n_high}  >90:{n_vhigh}  ({elapsed:.1f}s)")

    results.append({
        "name": name, "length": len(seq),
        "plddt_mean": round(mean_plddt, 1),
        "plddt_lt50": int(n_low), "plddt_50_70": int(n_mid),
        "plddt_70_90": int(n_high), "plddt_gt90": int(n_vhigh),
        "fold_time_s": round(elapsed, 1)
    })

# 保存结果
with open(ROOT / "work/round3/esmfold_results.json", "w") as f:
    json.dump(results, f, indent=2, ensure_ascii=False)

print(f"\n结果已保存到 work/round3/esmfold_results.json")

# 排序显示
print("\n" + "=" * 60)
print("pLDDT 排序 (mean, 常规标度 0-100)")
print("=" * 60)
results.sort(key=lambda x: x["plddt_mean"], reverse=True)
for r in results:
    bar = "█" * int(r["plddt_mean"] / 5)
    print(f"  {r['plddt_mean']:5.1f} {bar} {r['name']}")