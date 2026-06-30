"""
Round 5 P0-1: 扩展 LigandMPNN avGFP 设计 (多温度多种子)
==================================
策略:
  - 已确认 2WUR.pdb + chromo 区固定 = 100% 保留 TYG
  - 多温度 (T=0.1, 0.2, 0.3, 0.5) 多种子 (37, 137, 237) 生成更多多样性
  - 总目标: 100+ 高质量 LigandMPNN 候选

关键参数:
  - 固定 chromophore 区 (pos 30-42 in 2WUR, 含 TYG36-38) + 关键功能位
"""
import subprocess, json, sys, os, shutil
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
R5 = ROOT / "work" / "round5"
LMPNN = R5 / "LigandMPNN"
WORK = Path(r"C:\Temp\r5_lmpnn")
WORK.mkdir(parents=True, exist_ok=True)

shutil.copy(R5 / "pdbs" / "2WUR.pdb", WORK / "2WUR.pdb")

# 固定残基 (2WUR chromo 区 + 关键位)
FIXED_AV = " ".join([f"A{i}" for i in list(range(30, 43)) + [66, 118, 173, 192]])

# 多任务: 温度 + 种子组合
TASKS = []
for temp, n_batch in [("0.15", 5), ("0.2", 5), ("0.25", 5), ("0.4", 4)]:
    for seed in [137, 237]:
        TASKS.append({
            "name": f"av_lmpnn_T{temp.replace('.', '')}_s{seed}",
            "temp": temp,
            "seed": seed,
            "n_batch": n_batch,
        })

os.chdir(WORK)
print(f"工作目录: {os.getcwd()}\n共 {len(TASKS)} 任务\n")

for task in TASKS:
    out_dir = WORK / f"out_{task['name']}"
    out_dir.mkdir(exist_ok=True)
    print(f"\n[{task['name']}] T={task['temp']} seed={task['seed']} ({task['n_batch']*5} 序列)")
    
    cmd = [
        sys.executable, str(LMPNN / "run.py"),
        "--model_type", "ligand_mpnn",
        "--checkpoint_ligand_mpnn", str(LMPNN / "model_params" / "ligandmpnn_v_32_010_25.pt"),
        "--seed", str(task["seed"]),
        "--pdb_path", "2WUR.pdb",
        "--out_folder", str(out_dir),
        "--temperature", task["temp"],
        "--batch_size", "5",
        "--number_of_batches", str(task["n_batch"]),
        "--fixed_residues", FIXED_AV,
        "--save_stats", "1",
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600, cwd=str(WORK))
    if r.returncode != 0:
        print(f"  FAIL: {r.stderr[-200:]}")
    else:
        fa = out_dir / "seqs" / "2WUR.fa"
        if fa.exists():
            with open(fa) as f:
                content = f.read()
            blocks = [b for b in content.strip().split(">") if b.strip()]
            kept = 0
            for b in blocks[1:]:
                seq = b.split("\n", 1)[1].strip().replace("\n", "")
                if any(cb in seq for cb in ["TYG", "SYG", "GYG"]):
                    kept += 1
            print(f"  ✓ {len(blocks)-1} 序列, chromo 保留 {kept}/{len(blocks)-1}")

# 拷贝结果
proj_out = R5 / "lmpnn_expanded"
if proj_out.exists():
    shutil.rmtree(proj_out)
proj_out.mkdir()
for task in TASKS:
    src = WORK / f"out_{task['name']}"
    if src.exists():
        shutil.copytree(src, proj_out / task["name"])

print(f"\n✓ 结果拷贝到 {proj_out}")
