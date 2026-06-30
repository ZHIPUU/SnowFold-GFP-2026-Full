"""
Round 4 Step C v5: ProteinMPNN (终极版 - PDB放在工作目录根, 避免子目录嵌入name)
"""
import subprocess, json, sys, shutil, os
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
WORK = ROOT / "work" / "round4"

TMP_BASE = Path(r"C:\Temp\mpnn_work2")
TMP_BASE.mkdir(parents=True, exist_ok=True)

TMP_MPNN = TMP_BASE / "ProteinMPNN"
if not TMP_MPNN.exists():
    print("复制 ProteinMPNN 到临时目录...")
    shutil.copytree(WORK / "ProteinMPNN", TMP_MPNN)

# *** 关键: PDB 直接放在 TMP_BASE 根目录, 而不是子目录 ***
shutil.copy(WORK / "pdbs" / "2B3P.pdb", TMP_BASE / "2B3P.pdb")

helper = TMP_MPNN / "helper_scripts"
os.chdir(TMP_BASE)  # cwd = TMP_BASE

# Step 1: 用 --input_path "." 解析当前目录
print("Step 1: 解析 PDB...")
r = subprocess.run([
    sys.executable, str(helper / "parse_multiple_chains.py"),
    "--input_path", ".",
    "--output_path", "parsed.jsonl",
], capture_output=True, text=True, cwd=str(TMP_BASE))
print(" return:", r.returncode)
if r.returncode != 0:
    print(r.stderr); sys.exit(1)

# 验证 name
with open(TMP_BASE / "parsed.jsonl") as f:
    p = json.loads(f.readline())
print(f"  name = '{p['name']}'")
print(f"  seq len = {len(p.get('seq_chain_A',''))}")
seq_a = p.get('seq_chain_A', '')
for cb in ["TYG", "SYG", "GYG"]:
    if cb in seq_a:
        print(f"  {cb} at 1-based pos {seq_a.index(cb)+1}")
        break

# Step 2: 链分配
print("\nStep 2: 链分配...")
r = subprocess.run([
    sys.executable, str(helper / "assign_fixed_chains.py"),
    "--input_path", "parsed.jsonl",
    "--output_path", "chains.jsonl",
    "--chain_list", "A",
], capture_output=True, text=True, cwd=str(TMP_BASE))
print(" return:", r.returncode)

# Step 3: 固定位置
fixed_pos = "61 62 63 64 65 66 67 68 69 70 96 145 148 167 203 205 222"
print(f"\nStep 3: 固定位点 {fixed_pos}")
r = subprocess.run([
    sys.executable, str(helper / "make_fixed_positions_dict.py"),
    "--input_path", "parsed.jsonl",
    "--output_path", "fixed.jsonl",
    "--chain_list", "A",
    "--position_list", fixed_pos,
], capture_output=True, text=True, cwd=str(TMP_BASE))
print(" return:", r.returncode)

# Step 4: 跑 ProteinMPNN
configs = [("0.1", 30, "T01"), ("0.3", 30, "T03"), ("0.5", 20, "T05")]
for temp, n_seq, tag in configs:
    out_dir = f"out_{tag}"
    (TMP_BASE / out_dir).mkdir(exist_ok=True)
    print(f"\nStep 4: T={temp} ({n_seq} 序列)...")
    cmd = [
        sys.executable, str(TMP_MPNN / "protein_mpnn_run.py"),
        "--jsonl_path", "parsed.jsonl",
        "--chain_id_jsonl", "chains.jsonl",
        "--fixed_positions_jsonl", "fixed.jsonl",
        "--out_folder", out_dir,
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
        print(" STDERR tail:", r.stderr[-400:])
    else:
        fa_files = list((TMP_BASE / out_dir / "seqs").glob("*.fa"))
        for fa in fa_files:
            with open(fa) as f:
                content = f.read()
            n_seqs = sum(1 for l in content.split("\n") if l.startswith(">"))
            print(f"  ✓ {fa.name}: {n_seqs} 序列条目")

# 拷贝结果
proj_out = WORK / "mpnn_output_final"
if proj_out.exists():
    shutil.rmtree(proj_out)
proj_out.mkdir()
for tag in ["T01", "T03", "T05"]:
    src = TMP_BASE / f"out_{tag}"
    if src.exists():
        shutil.copytree(src, proj_out / tag)
print(f"\n✓ 结果拷贝到 {proj_out}")
