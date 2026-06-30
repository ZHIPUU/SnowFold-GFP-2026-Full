"""
R10: 迭代式 pLDDT 靶向修复管线
================================
基于 MPNN_T01_014 (最佳候选) 的残基级 pLDDT 分析：
 - 固定 pLDDT≥70 的高置信区域 + 生色团核心
 - 仅重设计低 pLDDT 区域
 - 多温度采样
 - ESMFold 验证
 - 评分排序
 - 迭代修复
"""
import json, torch, time, subprocess, sys, re
import numpy as np
from pathlib import Path

# ============ 配置 ============
ROOT = Path(r"D:\生信\2026Protein Design")
R10 = ROOT / "work" / "round10"
MPNN_DIR = Path(r"C:\proteinmpnn_r10")
MPNN_RUNNER = MPNN_DIR / "protein_mpnn_run.py"

# 最佳候选信息
BEST_NAME = "MPNN_T01_014"
BEST_SEQ = "MGKGDELFAGVVPVLVELDGDVNGHKFSVKGEGEGDASQGKLTLKFVCTTGELPVPWPTLVTTLTYGVQCFTRYPEHMKEHDFFKACMPEGYWRERTLKFKDDGTYKTRAEVKFEGDTLVNRIELKGTDFKEGGPILGHKIKYSYNSYNVYISPDKERNGIKATFTLRLDLEDGSTQLADVEELYTPIGEGPDELPRPHYLSVQNVLSKDPNEKRDHMVLEQFQSAGGIPAPMDELYK"
SEQ_LEN = len(BEST_SEQ)

# ============ ESMFold 推理函数 ============
def load_esmfold():
    from transformers import AutoTokenizer, EsmForProteinFolding
    t0 = time.time()
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
    model = EsmForProteinFolding.from_pretrained(
        "facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True).cuda()
    model.trunk.set_chunk_size(128)
    model.eval()
    print(f"  ESMFold 加载: {time.time()-t0:.1f}s")
    return tokenizer, model

def fold_sequence(seq, tokenizer, model):
    """返回 (mean_plddt, chromo_plddt, ptm, per_res_plddt)"""
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
    chromo_region = plddt_per_res[57:72]  # 0-based: 58-72
    cb_mean = float(chromo_region.mean())
    ptm = float(output.ptm.cpu().item())
    del tokens, output
    torch.cuda.empty_cache()
    return mean_plddt, cb_mean, ptm, plddt_per_res

def sort_score(ptm, plddt, chromo):
    """竞赛规则排序分 = pTM*0.55 + (pLDDT/100)*0.35 + (chromo/100)*0.10"""
    return ptm * 0.55 + (plddt / 100) * 0.35 + (chromo / 100) * 0.10

# ============ 步骤 1: 分析参考序列 ============
print("=" * 60)
print("R10 迭代修复管线 - 第 1 轮")
print("=" * 60)

print("\n[1/4] 加载 ESMFold 分析参考序列...")
tokenizer, model = load_esmfold()

t0 = time.time()
mean_plddt, chromo_plddt, ptm_val, per_res = fold_sequence(BEST_SEQ, tokenizer, model)
print(f"  参考: pLDDT={mean_plddt:.1f} chromo={chromo_plddt:.1f} pTM={ptm_val:.4f}")
print(f"  排序分: {sort_score(ptm_val, mean_plddt, chromo_plddt):.4f}")

# 确定固定/修复区域
absolute_core = {65, 66, 67, 96, 222}  # 生色团核心(1-based)
fixed_1based = set()
low_plddt_1based = []

for i in range(SEQ_LEN):
    res_idx = i + 1
    p = per_res[i]
    if p >= 70:
        fixed_1based.add(res_idx)
    elif p < 60:
        low_plddt_1based.append(res_idx)

# 确保核心位点也被固定
for r in absolute_core:
    fixed_1based.add(r)

fixed_str = " ".join(str(r) for r in sorted(fixed_1based))
print(f"  固定 {len(fixed_1based)} 个残基 (pLDDT≥70 + 核心)")
print(f"  需修复 {len(low_plddt_1based)} 个低 pLDDT 残基")
print(f"  固定列表: {fixed_str[:100]}...")

# ============ 步骤 2: ProteinMPNN 设计 ============
print("\n[2/4] 运行 ProteinMPNN 靶向修复...")

# 用 ESMFold 预测 PDB 作为输入
input_pdb = R10 / "pdbs" / "input.pdb"

# 使用预测 PDB + 固定残基列表
mpnn_out_dir = R10 / "mpnn_output"
mpnn_out_dir.mkdir(exist_ok=True)

# 多温度采样
temperatures = [0.1, 0.2, 0.3, 0.5]
seqs_per_temp = 40  # 每个温度 40 条 = 共 160 条

# 写固定残基 JSONL 文件 (ProteinMPNN 需要 JSONL 格式)
fixed_jsonl = mpnn_out_dir / "fixed_positions.jsonl"
with open(fixed_jsonl, "w") as f:
    f.write(json.dumps({"A": sorted(fixed_1based)}) + "\n")
print(f"  固定残基已写入: {fixed_jsonl}")

# 先解析 PDB 链
helper = MPNN_DIR / "helper_scripts" / "parse_multiple_chains.py"
r = subprocess.run([
    sys.executable, str(helper),
    "--input_path", str(R10 / "pdbs"),
    "--output_path", str(mpnn_out_dir / "parsed_chains.jsonl")
], capture_output=True, text=True)
print(f"  链解析: {'OK' if r.returncode == 0 else 'FAIL'}")
if r.returncode != 0:
    print(f"    {r.stderr[-300:]}")

all_candidates = []

for temp in temperatures:
    temp_out = mpnn_out_dir / f"T{str(temp).replace('.','')}"
    temp_out.mkdir(exist_ok=True)
    
    cmd = [
        sys.executable, str(MPNN_RUNNER),
        "--pdb_path", str(input_pdb),
        "--pdb_path_chains", "A",
        "--out_folder", str(temp_out),
        "--num_seq_per_target", str(seqs_per_temp),
        "--sampling_temp", str(temp),
        "--seed", str(int(time.time())),
        "--batch_size", "8",
        "--save_score", "1",
        "--model_name", "v_48_020",
        "--path_to_model_weights", str(MPNN_DIR / "vanilla_model_weights") + "/",
        "--fixed_positions_jsonl", str(fixed_jsonl),
    ]
    
    print(f"  T={temp}: 生成 {seqs_per_temp} 条...", end=" ", flush=True)
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=1200)
    
    if r.returncode != 0:
        print(f"FAIL ({time.time()-t0:.1f}s)")
        print(f"    {r.stderr[-500:]}")
        continue
    
    print(f"OK ({time.time()-t0:.1f}s)")
    
    # 读取生成的序列
    seq_file = temp_out / "seqs" / "input.fa"
    if seq_file.exists():
        content = seq_file.read_text()
        # 解析 FASTA
        seqs = []
        for block in content.strip().split(">"):
            if not block.strip():
                continue
            lines = block.strip().split("\n")
            header = lines[0]
            seq = "".join(lines[1:]).replace("\n", "").replace(" ", "")
            seqs.append(seq)
        
        for idx, seq in enumerate(seqs):
            all_candidates.append({
                "name": f"R10_T{str(temp).replace('.','')}_{idx+1:03d}",
                "temperature": temp,
                "seq": seq,
                "length": len(seq),
                "n_muts": sum(1 for a, b in zip(seq, BEST_SEQ) if a != b),
            })
        
        print(f"    读取 {len(seqs)} 条序列")

print(f"\n  总共生成 {len(all_candidates)} 条候选")

# ============ 步骤 3: ESMFold 验证 ============
print("\n[3/4] ESMFold 验证所有候选...")

for i, c in enumerate(all_candidates):
    seq = c["seq"]
    print(f"  [{i+1}/{len(all_candidates)}] {c['name']} ({c['n_muts']} muts)...", end=" ", flush=True)
    
    t0 = time.time()
    mp, cp, ptm_val, _ = fold_sequence(seq, tokenizer, model)
    c["plddt_mean"] = round(mp, 2)
    c["plddt_chromo"] = round(cp, 2)
    c["ptm"] = round(ptm_val, 4)
    c["sort_score"] = round(sort_score(ptm_val, mp, cp), 4)
    c["fold_time_s"] = round(time.time() - t0, 1)
    
    # 硬性门槛检查
    c["pass_ptm"] = ptm_val > 0.50
    c["pass_plddt"] = mp > 50.0
    c["pass_chromo"] = cp > 45.0
    c["all_pass"] = c["pass_ptm"] and c["pass_plddt"] and c["pass_chromo"]
    
    print(f"pLDDT={mp:.1f} chromo={cp:.1f} pTM={ptm_val:.4f} score={c['sort_score']:.4f} {'✅' if c['all_pass'] else '❌'}")

# ============ 步骤 4: 排序和选择 ============
print("\n[4/4] 排序和选择 Top 候选...")

# 按排序分降序排列
all_candidates.sort(key=lambda x: -x["sort_score"])

# 去重：相同的序列只保留第一个
seen_seqs = set()
unique_candidates = []
for c in all_candidates:
    if c["seq"] not in seen_seqs:
        seen_seqs.add(c["seq"])
        unique_candidates.append(c)

print(f"  去重后: {len(all_candidates)} → {len(unique_candidates)} 条唯一序列")

# 保存全部结果
with open(R10 / "r10_all_candidates.json", "w") as f:
    json.dump(unique_candidates, f, indent=2)

# Top 20
top20 = unique_candidates[:20]
with open(R10 / "r10_top20.json", "w") as f:
    json.dump(top20, f, indent=2)

print(f"\n  Top 20 已保存到 r10_top20.json")
print(f"\n{'='*60}")
print(f"Top 20 候选:")
print(f"{'排名':>4} {'名称':>20} {'突变':>4} {'pLDDT':>6} {'chromo':>6} {'pTM':>6} {'排序分':>7}")
print("-"*60)
for rank, c in enumerate(top20, 1):
    print(f"{rank:>4d} {c['name']:>20} {c['n_muts']:>4d} {c['plddt_mean']:>6.1f} {c['plddt_chromo']:>6.1f} {c['ptm']:>6.4f} {c['sort_score']:>7.4f}")

# 与参考序列对比
print(f"\n  参考 (MPNN_T01_014): pLDDT={mean_plddt:.1f} chromo={chromo_plddt:.1f} pTM={ptm_val:.4f} score={sort_score(ptm_val, mean_plddt, chromo_plddt):.4f}")

if top20:
    best = top20[0]
    delta_score = best["sort_score"] - sort_score(ptm_val, mean_plddt, chromo_plddt)
    print(f"  最佳提升: ΔpLDDT={best['plddt_mean']-mean_plddt:+.1f} Δchromo={best['plddt_chromo']-chromo_plddt:+.1f} ΔpTM={best['ptm']-ptm_val:+.4f} Δscore={delta_score:+.4f}")
    
    # 序列对比
    if best["seq"] != BEST_SEQ:
        print(f"\n  变化位点:")
        for i, (a, b) in enumerate(zip(best["seq"], BEST_SEQ)):
            if a != b:
                print(f"    Pos {i+1}: {b} → {a}")
