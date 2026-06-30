"""
Round 8: ProteinMPNN vanilla on Windows (CUDA + RTX 5080)
=========================================================
Model: v_48_020 (ProteinMPNN vanilla, Round 4 最佳模型)
固定: 生色团关键残基 (Y66, G67, R96, E222 + 核心骨架)
"""
import subprocess
import json
import sys
import shutil
import time
import os
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
R4 = ROOT / "work" / "round4"
R5 = ROOT / "work" / "round5"
R8 = ROOT / "work" / "round8"
R8.mkdir(parents=True, exist_ok=True)

PMPNN = Path(r"C:\LigandMPNN")
print(f"CUDA ready: {__import__('torch').cuda.is_available()}")

# 使用 Round 4 原始 PDB (固定残基格式)
pdb_dir = R4 / "pdbs"
assert pdb_dir.exists(), f"Missing: {pdb_dir}"

# 准备输出目录
OUT = R8 / "mpnn_designs"
OUT.mkdir(parents=True, exist_ok=True)

# 1. 准备固定残基: 生色团 + 催化 + 骨架核心 (共 17 个)
# 格式: A1 A2 A3 ... (chain_resid)
FIXED = "A58 A65 A66 A67 A68 A69 A70 A71 A72 A96 A148 A171 A203 A205 A206 A219 A222"
print(f"固定残基: {FIXED}")

# 2. 运行 ProteinMPNN (vanilla model, v_48_020)
# 多温度: T=0.1, 0.2, 0.3
TASKS = [
    ("T01", "0.1", 10),   # 保守
    ("T02", "0.2", 10),   # 中等
    ("T03", "0.3", 10),   # 多样
]

task_dirs = []

for tag, temp, n_batches in TASKS:
    tag_dir = Path(f"C:\\Temp\\r8_mpnn_{tag}")
    tag_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n--- ProteinMPNN T={temp} ---")
    t0 = time.time()
    
    # chdir 到 LigandMPNN 目录 (解决路径依赖)
    cwd = os.getcwd()
    os.chdir(str(PMPNN))
    
    # 复制 PDB 到工作目录 (LigandMPNN 需要)
    shutil.copy(str(pdb_dir / "2B3P.pdb"), "2B3P.pdb")
    
    # 中文路径编码问题，直接传 checkpoint 绝对路径
    ckpt = PMPNN / "model_params" / "proteinmpnn_v_48_020.pt"
    
    cmd = [
        sys.executable, "run.py",
        "--model_type", "protein_mpnn",
        "--checkpoint_protein_mpnn", str(ckpt),
        "--seed", "111",
        "--pdb_path", "2B3P.pdb",
        "--out_folder", str(tag_dir),
        "--temperature", temp,
        "--batch_size", "5",
        "--number_of_batches", str(n_batches),
        "--fixed_residues", FIXED,
        "--omit_AA", "C",
    ]
    
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    os.chdir(cwd)  # 恢复工作目录
    task_dirs.append(tag_dir)
    
    if r.returncode == 0:
        print(f"  ✅ {time.time()-t0:.0f}s")
    else:
        print(f"  ❌ rc={r.returncode}")
        print(f"  stderr: {r.stderr[:500]}")

# ============================================================
# 解析序列
# ============================================================
print(f"\n{'='*60}")
print("解析 ProteinMPNN 输出...")

WT = "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

candidates = []
seen = set()

for i, tag in enumerate([t[0] for t in TASKS]):
    tag_dir = task_dirs[i]
    fa_path = tag_dir / "seqs" / "2B3P.fa"
    if not fa_path.exists():
        print(f"  无输出: {tag}")
        # 尝试其他路径
        alt_path = list(tag_dir.rglob("*.fa"))
        if alt_path:
            fa_path = alt_path[0]
        else:
            continue
    
    with open(fa_path) as f:
        content = f.read()
    
    # ProteinMPNN 输出格式: >seq header\nSEQUENCE
    blocks = content.strip().split(">")
    for block in blocks:
        if not block.strip():
            continue
        lines = block.strip().split("\n", 1)
        if len(lines) < 2:
            continue
        seq = lines[1].strip().replace(" ", "").replace("\n", "")
        
        # 补全 X 为 WT 残基, 补 N 端 M
        seq_list = list(seq)
        for j, aa in enumerate(seq_list):
            if aa == "X":
                wt_idx = j + 1  # PDB offset
                if wt_idx < len(WT):
                    seq_list[j] = WT[wt_idx]
                else:
                    seq_list[j] = "A"
        full_seq = "M" + "".join(seq_list)
        
        # 合规验证
        if full_seq[0] != "M" or len(full_seq) < 220 or len(full_seq) > 250:
            continue
        if set(full_seq) - set("ACDEFGHIKLMNPQRSTVWY"):
            continue
        if full_seq not in seen:
            seen.add(full_seq)
            n_muts = sum(1 for a, b in zip(full_seq, WT) if a != b)
            candidates.append({
                "name": f"R8_PMPNN_{tag}_{len(candidates)+1:03d}",
                "scaffold": "sfGFP_MPNN",
                "n_muts": n_muts,
                "length": len(full_seq),
                "seq": full_seq,
                "notes": f"ProteinMPNN vanilla T={temp}"
            })

print(f"解析: {len(candidates)} 条 (去重)")
if candidates:
    # 多样性挑选
    candidates.sort(key=lambda x: x["n_muts"])
    print(f"  突变数范围: {candidates[0]['n_muts']}-{candidates[-1]['n_muts']}")

# 保存
with open(R8 / "candidates_mpnn_r8.json", "w", encoding="utf-8") as f:
    json.dump(candidates, f, indent=2, ensure_ascii=False)
print(f"\n保存: {R8 / 'candidates_mpnn_r8.json'}")
