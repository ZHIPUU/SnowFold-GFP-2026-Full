"""
Round 4 Step C v2: ProteinMPNN de novo (修正参数)
==================================
正确流程:
  1. parse_multiple_chains.py 生成 parsed_pdbs.jsonl
  2. assign_fixed_chains.py 指定要设计的链
  3. make_fixed_positions_dict.py 生成固定位置 jsonl
  4. protein_mpnn_run.py 运行设计
"""
import subprocess, json, sys
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
WORK = ROOT / "work" / "round4"
MPNN = WORK / "ProteinMPNN"
PDB_DIR = WORK / "pdbs_2b3p_only"  # 只放 sfGFP 一个
OUT_DIR = WORK / "mpnn_output"
OUT_DIR.mkdir(exist_ok=True)
PDB_DIR.mkdir(exist_ok=True)

# 复制 2B3P 到独立目录 (不让其他PDB干扰)
import shutil
shutil.copy(WORK / "pdbs" / "2B3P.pdb", PDB_DIR / "2B3P.pdb")

helper = MPNN / "helper_scripts"
parsed_jsonl = OUT_DIR / "parsed_pdbs.jsonl"
chains_jsonl = OUT_DIR / "assigned_chains.jsonl"
fixed_jsonl = OUT_DIR / "fixed_positions.jsonl"

# ============================================================
# Step 1: 解析 PDB
# ============================================================
print("Step 1: 解析 PDB...")
r = subprocess.run([
    sys.executable, str(helper / "parse_multiple_chains.py"),
    "--input_path", str(PDB_DIR),
    "--output_path", str(parsed_jsonl),
], capture_output=True, text=True)
if r.returncode != 0:
    print("FAIL:", r.stderr); sys.exit(1)
print(f"  生成 {parsed_jsonl}")
# 看一下解析结果
with open(parsed_jsonl) as f:
    parsed = [json.loads(l) for l in f]
for p in parsed:
    print(f"  {p['name']}: chains={list(p.keys())}, seq A len={len(p.get('seq_chain_A',''))}")

# ============================================================
# Step 2: 指定要设计的链 (A 链)
# ============================================================
print("\nStep 2: 指定设计链...")
r = subprocess.run([
    sys.executable, str(helper / "assign_fixed_chains.py"),
    "--input_path", str(parsed_jsonl),
    "--output_path", str(chains_jsonl),
    "--chain_list", "A",
], capture_output=True, text=True)
if r.returncode != 0:
    print("FAIL:", r.stderr); sys.exit(1)
print("  OK")

# ============================================================
# Step 3: 指定固定位置 (chromophore + 关键功能残基)
# ============================================================
# sfGFP chromophore TYG 在 pos 65-67 (1-based)
# 关键: R96 chromophore maturation, E222 proton wire, S205, H148 (proton acceptor)
# 也固定 chromophore 周围 8Å (近邻 ~10 个残基)
fixed_pos = "61 62 63 64 65 66 67 68 69 70 96 145 148 167 203 205 222"
print(f"\nStep 3: 固定位点 {fixed_pos}")
r = subprocess.run([
    sys.executable, str(helper / "make_fixed_positions_dict.py"),
    "--input_path", str(parsed_jsonl),
    "--output_path", str(fixed_jsonl),
    "--chain_list", "A",
    "--position_list", fixed_pos,
], capture_output=True, text=True)
if r.returncode != 0:
    print("FAIL:", r.stderr); sys.exit(1)
print("  OK")

# ============================================================
# Step 4: 运行 ProteinMPNN
# ============================================================
for temp, n_seq, tag in [("0.1", 30, "T01"), ("0.3", 30, "T03")]:
    out_sub = OUT_DIR / tag
    out_sub.mkdir(exist_ok=True)
    print(f"\nStep 4: ProteinMPNN T={temp} ({n_seq} 序列)...")
    cmd = [
        sys.executable, str(MPNN / "protein_mpnn_run.py"),
        "--jsonl_path", str(parsed_jsonl),
        "--chain_id_jsonl", str(chains_jsonl),
        "--fixed_positions_jsonl", str(fixed_jsonl),
        "--out_folder", str(out_sub),
        "--num_seq_per_target", str(n_seq),
        "--sampling_temp", temp,
        "--seed", "37",
        "--batch_size", "8",
        "--save_score", "1",
        "--model_name", "v_48_020",
        "--path_to_model_weights", str(MPNN / "vanilla_model_weights"),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        print("FAIL:", r.stderr[-1000:])
        print("STDOUT:", r.stdout[-500:])
    else:
        print(f"  OK (stdout tail: {r.stdout[-200:].strip()[:150]})")
        # 找到 .fa 文件
        fa_files = list((out_sub / "seqs").glob("*.fa"))
        for fa in fa_files:
            with open(fa) as f:
                content = f.read()
            n_lines = len([l for l in content.split("\n") if l.startswith(">")])
            print(f"  {fa.name}: {n_lines} 序列")

print("\n✓ 完成! 输出在 work/round4/mpnn_output/")
