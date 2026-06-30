"""
Round 4 Step C: ProteinMPNN de novo 序列生成
==================================
对 sfGFP (2B3P) 进行 inverse folding:
  - 固定 chromophore (T65, Y66, G67) 不变
  - 固定 β-barrel 关键骨架残基 (R96, E222)
  - 重设计其余位置
  - 温度 0.1 (保守) 和 0.5 (多样性)
  - 生成 50+ 候选
"""
import subprocess, json, sys, re
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
WORK = ROOT / "work" / "round4"
MPNN = WORK / "ProteinMPNN"
PDB_DIR = WORK / "pdbs"
OUT_DIR = WORK / "mpnn_output"
OUT_DIR.mkdir(exist_ok=True)

# 我们用 sfGFP 2B3P 作为模板 (现成的高质量结构)
target_pdb = PDB_DIR / "2B3P.pdb"

# sfGFP chromophore TYG 在 pos 65-67
# 关键功能位点: R96 (chromophore maturation), E222 (proton wire)
# 这些位点不重设计 (fixed)
fixed_positions = "65 66 67 96 222"  # ProteinMPNN 用 1-based

# ============================================================
# 准备 helper 命令
# ============================================================
helper_scripts = MPNN / "helper_scripts"

# 解析单链 PDB
parse_cmd = [
    sys.executable, str(helper_scripts / "parse_multiple_chains.py"),
    "--input_path", str(PDB_DIR),
    "--output_path", str(OUT_DIR / "parsed_chains.jsonl")
]
print("Step 1: 解析 PDB 链...")
print("  cmd:", " ".join(parse_cmd))
r = subprocess.run(parse_cmd, capture_output=True, text=True)
if r.returncode != 0:
    print("  STDERR:", r.stderr[-500:])
else:
    print("  OK")

# ============================================================
# 跑 ProteinMPNN (vanilla 模型 v_48_020, 中等噪声)
# ============================================================
# 关键参数:
#   --num_seq_per_target: 每个 target 生成多少序列 (50)
#   --sampling_temp: 采样温度 (用两个: 0.1 和 0.3)
#   --batch_size: 8
#   --pdb_path_chains: A (单链)
#   --pdb_path: sfGFP 2B3P
# ============================================================

run_cmds = [
    {
        "temp": "0.1",
        "out": OUT_DIR / "T01",
        "n_seq": 30,
    },
    {
        "temp": "0.3",
        "out": OUT_DIR / "T03",
        "n_seq": 30,
    },
]

for cfg in run_cmds:
    cfg["out"].mkdir(exist_ok=True)
    cmd = [
        sys.executable, str(MPNN / "protein_mpnn_run.py"),
        "--pdb_path", str(target_pdb),
        "--pdb_path_chains", "A",
        "--out_folder", str(cfg["out"]),
        "--num_seq_per_target", str(cfg["n_seq"]),
        "--sampling_temp", cfg["temp"],
        "--seed", "37",
        "--batch_size", "8",
        "--save_score", "1",
        "--model_name", "v_48_020",
        "--fixed_positions_list", fixed_positions,
        "--chain_id_jsonl", "",  # 默认A链
    ]
    print(f"\nStep 2: ProteinMPNN T={cfg['temp']} ({cfg['n_seq']} 序列)...")
    print("  cmd:", " ".join(cmd[:6]) + " ...")
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        print("  STDERR:", r.stderr[-1000:])
        print("  STDOUT:", r.stdout[-500:])
    else:
        print("  OK")
        print("  STDOUT tail:", r.stdout[-500:])
