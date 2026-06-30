"""
Round 4 Step 2: GPU ESMFold 批量评估
==================================
对所有 Round 4 候选 + Round 3 候选(对照) 做 GPU ESMFold 验证
速度: ~5s/序列 (GPU FP32), 比 CPU 快 ~26x
"""
import json, torch, time
from pathlib import Path
from transformers import AutoTokenizer, EsmForProteinFolding

ROOT = Path(r"D:\生信\2026Protein Design")
OUT = ROOT / "work" / "round4"

# 加载 Round 4 候选
with open(OUT / "candidates_round4.json", encoding="utf-8") as f:
    candidates = json.load(f)

print(f"Round 4 候选: {len(candidates)} 条")
for c in candidates:
    print(f"  [{c['role'][:20]:<20}] {c['name']:<32}  {c['n_muts']} mut")

# ============================================================
print("\n" + "=" * 70)
print(f"GPU: {torch.cuda.get_device_name(0)} ({torch.cuda.get_device_capability(0)})")
print(f"PyTorch: {torch.__version__}")
print(f"加载 ESMFold 模型...")
# ============================================================

tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained(
    "facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True
).cuda()
model.trunk.set_chunk_size(128)
model.eval()
print("✓ 模型就绪 (FP32, GPU)\n")

# ============================================================
# 批量推理
# ============================================================
results = []
for i, c in enumerate(candidates):
    name = c["name"]
    seq = c["seq"]
    print(f"[{i+1:>2}/{len(candidates)}] {name:<32} ({len(seq)} aa)...", end=" ", flush=True)
    t0 = time.time()

    with torch.no_grad():
        tokens = tokenizer([seq], return_tensors="pt", add_special_tokens=False)["input_ids"].cuda()
        output = model(tokens)

    # pLDDT 提取 (shape [1, L, 37]) 值域 [0,1] -> *100
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

    # 重点关注 chromophore 区域 (pos 64-67 ±5)
    cb_region = plddt_per_res[max(0,58):min(len(plddt_per_res),73)]
    cb_mean = float(cb_region.mean())

    # ptm (标量)
    try:
        ptm = float(output.ptm.cpu().item()) if hasattr(output, 'ptm') else None
    except Exception:
        ptm = None

    elapsed = time.time() - t0

    status_icon = "🟢" if mean_plddt >= 70 else ("🟡" if mean_plddt >= 50 else "🔴")
    print(f"pLDDT={mean_plddt:5.1f} cb={cb_mean:5.1f} {status_icon} ({elapsed:.1f}s)")

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

# ============================================================
# 保存 + 排序展示
# ============================================================
with open(OUT / "esmfold_round4.json", "w", encoding="utf-8") as f:
    json.dump(results, f, indent=2, ensure_ascii=False)

print(f"\n✓ 结果已保存到 work/round4/esmfold_round4.json")

# 排序按 pLDDT
print("\n" + "=" * 90)
print(f"{'name':<32} {'role':<22} {'mut':>3} {'pLDDT':>6} {'chromo_cb':>9} {'ptm':>6}")
print("=" * 90)
for r in sorted(results, key=lambda x: -x["plddt_mean"]):
    icon = "🟢" if r["plddt_mean"] >= 70 else ("🟡" if r["plddt_mean"] >= 50 else "🔴")
    ptm_str = f"{r['ptm']:.4f}" if r['ptm'] is not None else "  N/A"
    print(f"{r['name']:<32} {r['role'][:22]:<22} {r['n_muts']:>3} {r['plddt_mean']:>6.1f} {r['plddt_chromo_region']:>9.1f} {ptm_str:>6} {icon}")

# 按角色分组统计
print("\n按角色统计:")
from collections import defaultdict
by_role = defaultdict(list)
for r in results:
    by_role[r["role"]].append(r)
for role, items in by_role.items():
    items_sorted = sorted(items, key=lambda x: -x["plddt_mean"])
    print(f"\n  [{role}]:")
    for r in items_sorted:
        print(f"    pLDDT={r['plddt_mean']:5.1f}  {r['name']}")

print("\n下一步: 03_score_and_select.py 综合打分选 Top-6")
