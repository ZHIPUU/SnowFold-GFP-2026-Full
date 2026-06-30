"""
Round 4 Step C v4: ProteinMPNN (用相对路径修复 name 问题)
==================================
"""
import subprocess, json, sys, shutil, os
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
WORK = ROOT / "work" / "round4"
PROJECT_MPNN = WORK / "ProteinMPNN"

TMP_BASE = Path(r"C:\Temp\mpnn_work")
TMP_MPNN = TMP_BASE / "ProteinMPNN"
TMP_PDB = TMP_BASE / "pdbs"
TMP_OUT = TMP_BASE / "out"

shutil.copy(WORK / "pdbs" / "2B3P.pdb", TMP_PDB / "2B3P.pdb")

helper = TMP_MPNN / "helper_scripts"

# 关键: 切换工作目录到 TMP_BASE, 然后用相对路径
os.chdir(TMP_BASE)

# Step 1: 用相对路径 parse
print("Step 1: 解析 (用相对路径)...")
parsed = "out/parsed.jsonl"
r = subprocess.run([
    sys.executable, str(helper / "parse_multiple_chains.py"),
    "--input_path", "pdbs",
    "--output_path", parsed,
], capture_output=True, text=True, cwd=str(TMP_BASE))
print(" return:", r.returncode)
if r.returncode != 0:
    print(" STDERR:", r.stderr); sys.exit(1)

# 检查 name
with open(TMP_BASE / parsed) as f:
    p = json.loads(f.readline())
print(f"  解析的 name = '{p['name']}'")
print(f"  seq A len = {len(p.get('seq_chain_A',''))}")

# 找 chromophore
seq_a = p.get('seq_chain_A', '')
for cb in ["TYG", "SYG", "GYG"]:
    idx = seq_a.find(cb)
    if idx >= 0:
        print(f"  {cb} at 1-based pos {idx+1}-{idx+3}")
        break

# 关键位点 (基于实际 2B3P 序列)
# PDB 2B3P 可能从M开始, chromophore TYG 应在 65-67
fixed_positions = "61 62 63 64 65 66 67 68 69 70 96 145 148 167 203 205 222"
print(f"  固定位点: {fixed_positions}")

# Step 2: 链分配
chains_jsonl = "out/chains.jsonl"
r = subprocess.run([
    sys.executable, str(helper / "assign_fixed_chains.py"),
    "--input_path", parsed,
    "--output_path", chains_jsonl,
    "--chain_list", "A",
], capture_output=True, text=True, cwd=str(TMP_BASE))
print(f"Step 2 return: {r.returncode}")

# Step 3: 固定位置
fixed_jsonl = "out/fixed.jsonl"
r = subprocess.run([
    sys.executable, str(helper / "make_fixed_positions_dict.py"),
    "--input_path", parsed,
    "--output_path", fixed_jsonl,
    "--chain_list", "A",
    "--position_list", fixed_positions,
], capture_output=True, text=True, cwd=str(TMP_BASE))
print(f"Step 3 return: {r.returncode}")

# Step 4: ProteinMPNN
for temp, n_seq, tag in [("0.1", 30, "T01"), ("0.3", 30, "T03")]:
    out_sub = f"out/{tag}"
    (TMP_BASE / out_sub).mkdir(exist_ok=True)
    print(f"\nStep 4: T={temp} ({n_seq} 序列)...")
    cmd = [
        sys.executable, str(TMP_MPNN / "protein_mpnn_run.py"),
        "--jsonl_path", parsed,
        "--chain_id_jsonl", chains_jsonl,
        "--fixed_positions_jsonl", fixed_jsonl,
        "--out_folder", out_sub,
        "--num_seq_per_target", str(n_seq),
        "--sampling_temp", temp,
        "--seed", "37",
        "--batch_size", "8",
        "--save_score", "1",
        "--model_name", "v_48_020",
        "--path_to_model_weights", str(TMP_MPNN / "vanilla_model_weights"),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600, cwd=str(TMP_BASE))
    print(f"  return: {r.returncode}")
    if r.returncode != 0:
        print(" STDERR tail:", r.stderr[-600:])
    else:
        fa_files = list((TMP_BASE / out_sub / "seqs").glob("*.fa"))
        for fa in fa_files:
            with open(fa) as f:
                content = f.read()
            n_seqs = sum(1 for l in content.split("\n") if l.startswith(">"))
            print(f"  {fa.name}: {n_seqs} 序列")

# 拷贝结果回项目
proj_out = WORK / "mpnn_output_v3"
if proj_out.exists():
    shutil.rmtree(proj_out)
shutil.copytree(TMP_BASE / "out", proj_out)
print(f"\n✓ 拷贝到 {proj_out}")
