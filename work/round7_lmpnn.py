"""
Round 7 Phase 2: LigandMPNN chromophore-aware design
=====================================================
核心: LigandMPNN 显式建模 CRO chromophore 作为配体
比 ProteinMPNN (无视 chromophore) 更适合 GFP

固定策略 (最少固定):
  - CRO: chromophore 自身 (LigandMPNN 自动识别 HETATM)
  - Only Y66: 生色团三联体 (固定), R96 (成熟催化), E222 (质子传递)
  - 其余全部开放给 LigandMPNN 设计

输出: 50+ 候选 → ESMFold round-trip → 严格过滤
"""
import subprocess, json, sys, shutil, time, os
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
R5 = ROOT / "work" / "round5"
LMPNN = R5 / "LigandMPNN"
OUT = ROOT / "work" / "round7"
OUT.mkdir(parents=True, exist_ok=True)

WORK = Path(r"C:\Temp\r7_lmpnn")
if WORK.exists():
    shutil.rmtree(WORK)
WORK.mkdir(parents=True, exist_ok=True)

# 复制 PDB
shutil.copy(R5 / "pdbs" / "2B3P.pdb", WORK / "2B3P.pdb")

print("=" * 80)
print("Round 7: LigandMPNN - chromophore-aware GFP design")
print("=" * 80)

# 检查权重
weights_dir = LMPNN / "model_params"
ckpt = weights_dir / "ligandmpnn_v_32_020_25.pt"
assert ckpt.exists(), f"权重缺失: {ckpt}"
print(f"权重: {ckpt.name} ({ckpt.stat().st_size/1e6:.0f}MB)")

# ============================================================
# 任务: 多温度 + 多批次
# ============================================================
# 固定仅关键位 (LigandMPNN 格式: chain_residue)
# A66 = Y66 (chromophore 中间), A96 = R96, A222 = E222
FIXED = "A66 A96 A222"

TASKS = [
    ("T01", "0.1", 10, 5),   # T=0.1, 10 batches × 5 = 50 seqs
    ("T03", "0.3", 10, 5),   # T=0.3
    ("T05", "0.5", 10, 5),   # T=0.5
]

for tag, temp, n_batches, batch_size in TASKS:
    out_dir = WORK / tag
    out_dir.mkdir(exist_ok=True)
    
    print(f"\n--- LigandMPNN T={temp} ({n_batches}×{batch_size}={n_batches*batch_size} seqs) ---")
    t0 = time.time()
    
    cmd = [
        sys.executable, str(LMPNN / "run.py"),
        "--model_type", "ligand_mpnn",
        "--seed", "111",
        "--pdb_path", "2B3P.pdb",
        "--out_folder", str(out_dir),
        "--temperature", temp,
        "--batch_size", str(batch_size),
        "--number_of_batches", str(n_batches),
        "--fixed_residues", FIXED,
        "--save_stats", "1",
        "--checkpoint_ligand_mpnn", str(ckpt),
    ]
    
    # 切换到工作目录 (LigandMPNN 相对路径依赖)
    cwd = os.getcwd()
    os.chdir(str(WORK))
    
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    os.chdir(cwd)  # 恢复目录
    elapsed = time.time() - t0
    if r.returncode == 0:
        print(f"  ✅ {elapsed:.0f}s")
    else:
        print(f"  ❌ rc={r.returncode}")
        print(f"  stderr: {r.stderr[:500]}")

# ============================================================
# 解析输出
# ============================================================
print(f"\n{'='*80}")
print("解析 LigandMPNN 输出...")

LMPNN_WT = "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

candidates = []
seen = set()

for tag, temp, _, _ in TASKS:
    fa_path = WORK / tag / "seqs" / "2B3P.fa"
    if not fa_path.exists():
        print(f"  无输出: {tag}")
        continue
    
    with open(fa_path) as f:
        content = f.read()
    
    # LigandMPNN 输出中固定残基也是 X (从 PDB 读取)
    # 需要用 WT 序列的对应位置填回
    blocks = content.strip().split(">")
    for i, block in enumerate(blocks):
        if not block.strip() or i == 0:  # skip header/empty
            continue
        lines = block.strip().split("\n", 1)
        if len(lines) < 2:
            continue
        seq = lines[1].strip().replace("\n", "")
        
        # LigandMPNN 输出从 PDB 序列 (231 aa, 从 S 开始)
        # 需要对齐到 sfGFP WT (238 aa, 从 M 开始)
        # PDB pos 1 = WT pos 2 = 'S'
        # 补 N 端 M, 填 X
        seq_list = list(seq)
        for j, aa in enumerate(seq_list):
            if aa == "X":
                wt_idx = j + 1  # PDB offset: PDB j = WT j+1
                if wt_idx < len(LMPNN_WT):
                    seq_list[j] = LMPNN_WT[wt_idx]
        full_seq = "M" + "".join(seq_list)
        
        if full_seq[0] != "M" or len(full_seq) < 220 or len(full_seq) > 250:
            continue
        if set(full_seq) - set("ACDEFGHIKLMNPQRSTVWY"):
            continue
        if full_seq not in seen:
            seen.add(full_seq)
            n_muts = sum(1 for a, b in zip(full_seq, LMPNN_WT) if a != b)
            candidates.append({
                "name": f"R7_LMPNN_{tag}_{i:03d}",
                "scaffold": "sfGFP_LMPNN",
                "n_muts": n_muts,
                "length": len(full_seq),
                "seq": full_seq,
                "notes": f"LigandMPNN T={temp} fixed={FIXED}"
            })

print(f"解析完成: {len(candidates)} 条候选 (去重)")

# 保存
with open(OUT / "candidates_lmpnn_r7.json", "w", encoding="utf-8") as f:
    json.dump(candidates, f, indent=2, ensure_ascii=False)

# 统计
muts = [c["n_muts"] for c in candidates]
if muts:
    print(f"  突变数: avg={sum(muts)/len(muts):.0f} range={min(muts)}-{max(muts)}")
    candidates.sort(key=lambda x: x["n_muts"])
    print(f"  最少突变 Top-5:")
    for c in candidates[:5]:
        print(f"    {c['name']}: {c['n_muts']} muts, len={c['length']}")
else:
    print("  ❌ 无候选生成!")
    print(f"  检查 LigandMPNN 输出路径: {WORK}")
