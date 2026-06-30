"""ESMFold pLDDT 验证所有 Round 3 候选"""
import json, torch
import numpy as np
from pathlib import Path
from transformers import AutoTokenizer, EsmForProteinFolding

ROOT = Path(r"D:\生信\2026Protein Design")

# 加载候选
with open(ROOT / "work/round3/candidates_round3.json") as f:
    candidates = json.load(f)

print(f"加载 ESMFold 模型...")
tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True)
model = model.cuda()
# Use FP32 for accurate folding (FP16 degrades pLDDT)
model.trunk.set_chunk_size(128)  # balance memory/quality
print("模型就绪\n")

results = []
for i, c in enumerate(candidates):
    name = c["name"]
    seq = c["seq"]
    print(f"[{i+1}/{len(candidates)}] Folding {name} ({len(seq)} aa)...")
    
    with torch.no_grad():
        tokens = tokenizer([seq], return_tensors="pt", add_special_tokens=False)["input_ids"].cuda()
        output = model(tokens)
    
    # 提取 pLDDT (shape: [1, L, 37]) 值域 [0,1]
    plddt_raw = output["plddt"].cpu().numpy()[0]          # (L, 37)
    atom_mask = output["atom37_atom_exists"].cpu().numpy()[0]  # (L, 37) bool
    
    # 转成常规 pLDDT [0,100] 并用 mask 排除虚拟原子
    plddt_scaled = plddt_raw * 100.0                       # (L, 37)
    # 对每个残基只统计实际存在的原子
    masked_sums = (plddt_scaled * atom_mask).sum(axis=1)   # (L,)
    masked_counts = atom_mask.sum(axis=1).astype(float)    # (L,)
    masked_counts[masked_counts == 0] = 1                  # 防除零
    plddt_per_res = masked_sums / masked_counts            # (L,) per-residue pLDDT [0,100]
    mean_plddt = float(plddt_per_res.mean())
    
    # 计算 pLDDT 分布 (常规标度)
    n_low = int((plddt_per_res < 50).sum().item())
    n_mid = int(((plddt_per_res >= 50) & (plddt_per_res < 70)).sum().item())
    n_high = int(((plddt_per_res >= 70) & (plddt_per_res < 90)).sum().item())
    n_vhigh = int((plddt_per_res >= 90).sum().item())
    
    status = "✓" if mean_plddt >= 70 else "✗"
    print(f"  pLDDT mean={mean_plddt:.1f} {status}  (<50:{n_low}  50-70:{n_mid}  70-90:{n_high}  >90:{n_vhigh})")
    
    results.append({
        "name": name, "length": len(seq),
        "plddt_mean": round(mean_plddt, 1),
        "plddt_lt50": int(n_low), "plddt_50_70": int(n_mid),
        "plddt_70_90": int(n_high), "plddt_gt90": int(n_vhigh)
    })
    
    # 清理显存
    del tokens, output
    torch.cuda.empty_cache()

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
