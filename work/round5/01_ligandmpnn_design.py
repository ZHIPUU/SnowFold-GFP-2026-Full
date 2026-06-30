"""
Round 5 P0-1: LigandMPNN GFP 设计 (chromophore-aware)
==================================================
关键创新:
  - 2B3P/2WUR PDB 含 CRO chromophore (HETATM)
  - LigandMPNN 把 CRO 当配体, 自动给周围残基提供化学约束
  - 比 Round 4 ProteinMPNN (无视 chromophore) 更智能

固定残基策略:
  - chromophore CRO 本身 (LigandMPNN 自动识别)
  - 关键功能位 R96 (chromophore maturation), E222 (proton wire), H148, T203
  - 8Å 内的近邻残基 (LigandMPNN 自动加权)
"""
import subprocess, json, sys, os, shutil
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
R5 = ROOT / "work" / "round5"
LMPNN = R5 / "LigandMPNN"

# 用纯英文工作目录避免中文路径 bug (Round 4 教训)
WORK = Path(r"C:\Temp\r5_lmpnn")
WORK.mkdir(parents=True, exist_ok=True)

# 复制 PDB
shutil.copy(R5 / "pdbs" / "2B3P.pdb", WORK / "2B3P.pdb")
shutil.copy(R5 / "pdbs" / "2WUR.pdb", WORK / "2WUR.pdb")

# 模型权重检查
weights_dir = LMPNN / "model_params"
required = {
    "ligand_mpnn": weights_dir / "ligandmpnn_v_32_010_25.pt",  # 中等噪声
    "soluble_mpnn": weights_dir / "solublempnn_v_48_020.pt",
}

print("模型权重检查:")
for name, p in required.items():
    if p.exists() and p.stat().st_size > 5_000_000:
        print(f"  ✓ {name}: {p.stat().st_size/1024/1024:.2f}MB")
    else:
        print(f"  ✗ {name}: MISSING or incomplete")

# ============================================================
# 任务清单
# ============================================================
# 对 2B3P sfGFP: chromophore CRO 在 chain A pos 66 (含 SER65/TYR66/GLY67 三联体)
# 关键功能位 (要固定): 96, 145, 148, 167, 203, 205, 222
# fixed_residues 用 LigandMPNN 格式: "A96 A148 A203 A222"

TASKS = [
    {
        "name": "sfGFP_lmpnn_T01",
        "pdb": "2B3P.pdb",
        "model": "ligand_mpnn",
        "checkpoint": "model_params/ligandmpnn_v_32_010_25.pt",
        "temperature": "0.1",
        "n_seq_per_batch": 5,
        "n_batches": 6,  # 30 总
        "fixed": "A66 A96 A148 A203 A205 A222",
        "out": "out_sf_lmpnn_T01",
    },
    {
        "name": "sfGFP_lmpnn_T03",
        "pdb": "2B3P.pdb",
        "model": "ligand_mpnn",
        "checkpoint": "model_params/ligandmpnn_v_32_010_25.pt",
        "temperature": "0.3",
        "n_seq_per_batch": 5,
        "n_batches": 4,  # 20 总
        "fixed": "A66 A96 A148 A203 A205 A222",
        "out": "out_sf_lmpnn_T03",
    },
    {
        "name": "avGFP_lmpnn_T01",
        "pdb": "2WUR.pdb",
        "model": "ligand_mpnn",
        "checkpoint": "model_params/ligandmpnn_v_32_010_25.pt",
        "temperature": "0.1",
        "n_seq_per_batch": 5,
        "n_batches": 6,
        "fixed": "A66 A96 A148 A203 A205 A222",  # avGFP 编号可能略偏, 后续可调
        "out": "out_av_lmpnn_T01",
    },
    # SolubleMPNN 对照 (无 ligand context)
    {
        "name": "sfGFP_soluble_T01",
        "pdb": "2B3P.pdb",
        "model": "soluble_mpnn",
        "checkpoint": "model_params/solublempnn_v_48_020.pt",
        "temperature": "0.1",
        "n_seq_per_batch": 5,
        "n_batches": 4,  # 20 总
        "fixed": "A66 A96 A148 A203 A205 A222",
        "out": "out_sf_soluble_T01",
    },
]

# ============================================================
# 准备脚本到 WORK 目录
# ============================================================
# 复制 LigandMPNN 项目到工作区 (用 link 而非 copy 节省空间)
LMPNN_LINK = WORK / "LigandMPNN"
if not LMPNN_LINK.exists():
    try:
        os.symlink(str(LMPNN), str(LMPNN_LINK), target_is_directory=True)
        print(f"\n创建符号链接 {LMPNN_LINK}")
    except OSError:
        # Windows 没管理员权限可能失败, 改用 mklink junction
        subprocess.run(["cmd", "/c", "mklink", "/J", str(LMPNN_LINK), str(LMPNN)], capture_output=True)
        print(f"\n创建 junction {LMPNN_LINK}")

os.chdir(WORK)
print(f"\n工作目录: {os.getcwd()}\n")

# ============================================================
# 跑每个任务
# ============================================================
for task in TASKS:
    out_dir = WORK / task["out"]
    out_dir.mkdir(exist_ok=True)
    
    print(f"\n{'='*70}")
    print(f"任务: {task['name']}  ({task['n_seq_per_batch'] * task['n_batches']} 序列)")
    print('='*70)
    
    # 仅 ligand_mpnn 模式需要的检查点参数
    ckpt_args = []
    if task["model"] == "ligand_mpnn":
        ckpt_args = ["--checkpoint_ligand_mpnn", str(LMPNN / task["checkpoint"])]
    elif task["model"] == "soluble_mpnn":
        ckpt_args = ["--checkpoint_soluble_mpnn", str(LMPNN / task["checkpoint"])]
    
    # 检查 checkpoint 是否存在
    if ckpt_args:
        ckpt_path = Path(ckpt_args[1])
        if not ckpt_path.exists():
            print(f"  ⚠ checkpoint 不存在: {ckpt_path}, 跳过此任务")
            continue
    
    cmd = [
        sys.executable, str(LMPNN / "run.py"),
        "--model_type", task["model"],
        "--seed", "37",
        "--pdb_path", task["pdb"],
        "--out_folder", str(out_dir),
        "--temperature", task["temperature"],
        "--batch_size", str(task["n_seq_per_batch"]),
        "--number_of_batches", str(task["n_batches"]),
        "--fixed_residues", task["fixed"],
        "--save_stats", "1",
    ] + ckpt_args
    
    print(f"  cmd: {cmd[0]} {cmd[1]} --model_type {task['model']} --temperature {task['temperature']}")
    print(f"  固定残基: {task['fixed']}")
    
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=900, cwd=str(WORK))
    print(f"  return: {r.returncode}")
    if r.returncode != 0:
        print(f"  STDERR (tail):\n{r.stderr[-800:]}")
    else:
        # 检查输出
        seqs_dir = out_dir / "seqs"
        if seqs_dir.exists():
            fa_files = list(seqs_dir.glob("*.fa"))
            for fa in fa_files:
                with open(fa) as f:
                    content = f.read()
                n = sum(1 for l in content.split("\n") if l.startswith(">"))
                print(f"  ✓ {fa.name}: {n} 序列")

# ============================================================
# 拷贝结果回项目目录
# ============================================================
out_proj = R5 / "lmpnn_output"
if out_proj.exists():
    shutil.rmtree(out_proj)
out_proj.mkdir()
for task in TASKS:
    src = WORK / task["out"]
    if src.exists():
        shutil.copytree(src, out_proj / task["name"])

print(f"\n\n✓ 结果拷贝到 {out_proj}")
