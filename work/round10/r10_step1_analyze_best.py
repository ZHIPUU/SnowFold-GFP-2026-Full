"""
R10 Step 1: 分析最佳候选 MPNN_T01_014 的残基级 pLDDT
=================================================
1. ESMFold 预测结构 + 获取每个残基的 pLDDT
2. 保存预测 PDB（给 ProteinMPNN 用）
3. 识别低 pLDDT 区域 (<60)
4. 输出固位置列表
"""
import json, torch, time, numpy as np
from pathlib import Path

# ============ 配置 ============
ROOT = Path(r"D:\生信\2026Protein Design")
R10 = ROOT / "work" / "round10"
R10.mkdir(exist_ok=True)

# MPNN_T01_014 序列（从 final_top6.json 提取）
BEST_SEQ = "MGKGDELFAGVVPVLVELDGDVNGHKFSVKGEGEGDASQGKLTLKFVCTTGELPVPWPTLVTTLTYGVQCFTRYPEHMKEHDFFKACMPEGYWRERTLKFKDDGTYKTRAEVKFEGDTLVNRIELKGTDFKEGGPILGHKIKYSYNSYNVYISPDKERNGIKATFTLRLDLEDGSTQLADVEELYTPIGEGPDELPRPHYLSVQNVLSKDPNEKRDHMVLEQFQSAGGIPAPMDELYK"
SEQ_LEN = len(BEST_SEQ)
print(f"MPNN_T01_014: {SEQ_LEN} aa")

# ============ ESMFold 预测 ============
from transformers import AutoTokenizer, EsmForProteinFolding

print("加载 ESMFold 模型...")
t0 = time.time()
tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained(
    "facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True).cuda()
model.trunk.set_chunk_size(128)
model.eval()
print(f"  模型加载完成 ({time.time()-t0:.1f}s)")

print("运行 ESMFold 预测...")
t0 = time.time()
with torch.no_grad():
    tokens = tokenizer([BEST_SEQ], return_tensors="pt", add_special_tokens=False)["input_ids"].cuda()
    output = model(tokens)
elapsed = time.time() - t0
print(f"  预测完成 ({elapsed:.1f}s)")

# ============ 提取 pLDDT ============
plddt_raw = output.plddt.cpu().numpy()[0]       # [L, 37_atoms]
atom_mask = output.atom37_atom_exists.cpu().numpy()[0]  # [L, 37]
plddt_scaled = plddt_raw * 100.0
masked_sums = (plddt_scaled * atom_mask).sum(axis=1)
masked_counts = atom_mask.sum(axis=1).astype(float)
masked_counts[masked_counts == 0] = 1
plddt_per_res = masked_sums / masked_counts     # [L] 每个残基的 pLDDT

mean_plddt = float(plddt_per_res.mean())
chromo_region = plddt_per_res[57:72]  # 0-based: 残基 58-72
cb_mean = float(chromo_region.mean())
ptm = float(output.ptm.cpu().item())

print(f"\n=== 总体指标 ===")
print(f"  全局 pLDDT: {mean_plddt:.1f}")
print(f"  生色团(58-72) pLDDT: {cb_mean:.1f}")
print(f"  pTM: {ptm:.4f}")

# ============ 保存预测 PDB ============
pdb_path = R10 / "mpnn_t01_014_pred.pdb"
with open(pdb_path, "w") as f:
    f.write(model.output_to_pdb(output)[0])
print(f"\n预测 PDB 已保存: {pdb_path}")

# ============ pLDDT 分析 ============
print(f"\n=== 每个残基 pLDDT 分析 ===")
print(f"{'残基':>5} {'AA':>3} {'pLDDT':>6} {'状态':>8}")
print("-"*30)

low_plddt_residues = []   # pLDDT < 60 需要修复
medium_plddt_residues = []  # 60 ≤ pLDDT < 70 
fixed_residues = []        # pLDDT ≥ 70 + 生色团 + 核心功能位点

# 绝对核心：生色团三联体 + 成熟关键位点
absolute_core = {65, 66, 67, 96, 222}  # TYG, R96, E222 (1-based)

for i in range(SEQ_LEN):
    res_idx_1based = i + 1
    aa = BEST_SEQ[i]
    p = plddt_per_res[i]
    
    if p >= 70:
        status = "固定✅"
        fixed_residues.append(res_idx_1based)
    elif p >= 60:
        status = "边缘⚠️"
        medium_plddt_residues.append(res_idx_1based)
    else:
        status = "低❌"
        low_plddt_residues.append(res_idx_1based)

    print(f"{res_idx_1based:>5d} {aa:>3} {p:>6.1f} {status:>8}")

# 确保绝对核心位点在 fixed 列表
for core_res in absolute_core:
    if core_res not in fixed_residues:
        fixed_residues.append(core_res)

fixed_residues = sorted(set(fixed_residues))

print(f"\n=== 统计 ===")
print(f"  固定(pLDDT≥70+核心): {len(fixed_residues)} 个残基")
print(f"  边缘(60-70): {len(medium_plddt_residues)} 个残基")
print(f"  低pLDDT(<60,需修复): {len(low_plddt_residues)} 个残基")

print(f"\n=== 低 pLDDT 区域 ===")
if low_plddt_residues:
    # 找连续区域
    regions = []
    start = low_plddt_residues[0]
    end = low_plddt_residues[0]
    for r in low_plddt_residues[1:]:
        if r == end + 1:
            end = r
        else:
            regions.append((start, end))
            start = r
            end = r
    regions.append((start, end))
    for s, e in regions:
        span = f"{s}-{e}" if s != e else str(s)
        # 打印这段的序列
        seq_span = BEST_SEQ[s-1:e]
        print(f"  残基 {span}: {seq_span}")

print(f"\n=== 固定残基列表 (ProteinMPNN 1-based) ===")
print(" ".join(str(r) for r in fixed_residues))

# ============ 保存分析结果 ============
result = {
    "name": "MPNN_T01_014",
    "seq": BEST_SEQ,
    "length": SEQ_LEN,
    "mean_plddt": round(mean_plddt, 2),
    "chromo_plddt": round(cb_mean, 2),
    "ptm": round(ptm, 4),
    "per_residue_plddt": [round(float(p), 2) for p in plddt_per_res],
    "fixed_residues_1based": fixed_residues,
    "low_plddt_residues_1based": low_plddt_residues,
    "medium_plddt_residues_1based": medium_plddt_residues,
}

with open(R10 / "t01_014_analysis.json", "w") as f:
    json.dump(result, f, indent=2)
print(f"\n分析结果已保存: {R10 / 't01_014_analysis.json'}")
