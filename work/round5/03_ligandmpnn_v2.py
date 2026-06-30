"""
Round 5 P0-1 修正版: LigandMPNN 设计 — 正确固定 chromophore 区
==================================
关键修正 (vs v1):
  - 2B3P 中 CRO 占 pos 66 (融合 SER65-TYR66-GLY67)
  - 固定 pos 60-72 整个 chromophore 区域 (包括邻近 helix)
  - 还固定 R96, E222, H148, T203 等关键功能位

新设计:
  - 减小 fixed_residues 范围, 只保留 chromophore 周围必要位
  - 允许 LigandMPNN 设计 N/C 端表面 + β-barrel
"""
import subprocess, json, sys, os, shutil
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
R5 = ROOT / "work" / "round5"
LMPNN = R5 / "LigandMPNN"

WORK = Path(r"C:\Temp\r5_lmpnn")
WORK.mkdir(parents=True, exist_ok=True)

# 复制 PDB
for pdb in ["2B3P.pdb", "2WUR.pdb"]:
    shutil.copy(R5 / "pdbs" / pdb, WORK / pdb)

# 固定残基策略 (basic_2B3P)
# - PDB 起始残基 pos 2 (缺 M), chromophore CRO 在 PDB pos 66
# - chromophore 周围 ±6 = pos 60-72
# - 关键功能位: R96, H148, T203, E222
# - 也固定 N端 M-S 区 (避免 ProteinMPNN 修改 N 端导致问题)
# 注意 LigandMPNN 用 "A66" 不是 "66"
FIXED_SF = " ".join([f"A{i}" for i in list(range(60, 73)) + [96, 148, 203, 222]])
# avGFP 2WUR chromophore TYG 在 1-based pos 36 (因为起始残基 pos 1 = T)
# 实际 avGFP chromo pos 65-67, 但 2WUR PDB 偏移导致 pos 36
# 同样固定 chromophore ± 6 + 关键位 (偏移 30)
FIXED_AV = " ".join([f"A{i}" for i in list(range(30, 43)) + [66, 118, 173, 192]])

TASKS = [
    {
        "name": "sfGFP_lmpnn_v2_T01",
        "pdb": "2B3P.pdb",
        "ckpt": "ligandmpnn_v_32_010_25.pt",
        "temp": "0.1",
        "n_batches": 6,
        "fixed": FIXED_SF,
        "out": "out_sf_lmpnn_v2_T01",
    },
    {
        "name": "sfGFP_lmpnn_v2_T03",
        "pdb": "2B3P.pdb",
        "ckpt": "ligandmpnn_v_32_010_25.pt",
        "temp": "0.3",
        "n_batches": 4,
        "fixed": FIXED_SF,
        "out": "out_sf_lmpnn_v2_T03",
    },
    {
        "name": "avGFP_lmpnn_v2_T01",
        "pdb": "2WUR.pdb",
        "ckpt": "ligandmpnn_v_32_010_25.pt",
        "temp": "0.1",
        "n_batches": 6,
        "fixed": FIXED_AV,
        "out": "out_av_lmpnn_v2_T01",
    },
]

os.chdir(WORK)
print(f"工作目录: {os.getcwd()}\n")

for task in TASKS:
    out_dir = WORK / task["out"]
    out_dir.mkdir(exist_ok=True)
    print(f"\n{'='*70}")
    print(f"任务: {task['name']}  ({task['n_batches'] * 5} 序列)")
    print('='*70)
    print(f"  固定残基: {task['fixed'][:100]}...")
    
    cmd = [
        sys.executable, str(LMPNN / "run.py"),
        "--model_type", "ligand_mpnn",
        "--checkpoint_ligand_mpnn", str(LMPNN / "model_params" / task["ckpt"]),
        "--seed", "37",
        "--pdb_path", task["pdb"],
        "--out_folder", str(out_dir),
        "--temperature", task["temp"],
        "--batch_size", "5",
        "--number_of_batches", str(task["n_batches"]),
        "--fixed_residues", task["fixed"],
        "--save_stats", "1",
    ]
    
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=900, cwd=str(WORK))
    print(f"  return: {r.returncode}")
    if r.returncode != 0:
        print(f"  STDERR: {r.stderr[-500:]}")
    else:
        seqs_dir = out_dir / "seqs"
        if seqs_dir.exists():
            for fa in seqs_dir.glob("*.fa"):
                with open(fa) as f:
                    content = f.read()
                n = sum(1 for l in content.split("\n") if l.startswith(">"))
                # 验证 chromophore 保留
                blocks = [b for b in content.strip().split(">") if b.strip()]
                kept = 0
                for b in blocks[1:]:
                    seq = b.split("\n", 1)[1].strip().replace("\n","")
                    if any(cb in seq for cb in ["TYG","SYG","GYG"]):
                        kept += 1
                print(f"  ✓ {fa.name}: {n} 序列, chromophore 保留 {kept}/{n-1}")

# 拷贝结果
proj_out = R5 / "lmpnn_output_v2"
if proj_out.exists():
    shutil.rmtree(proj_out)
proj_out.mkdir()
for task in TASKS:
    src = WORK / task["out"]
    if src.exists():
        shutil.copytree(src, proj_out / task["name"])
print(f"\n✓ 结果拷贝到 {proj_out}")
