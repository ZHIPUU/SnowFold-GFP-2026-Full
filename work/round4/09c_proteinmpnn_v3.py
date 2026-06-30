"""
Round 4 Step C v3: ProteinMPNN de novo (避开中文路径)
==================================
方案: 把 ProteinMPNN 工作目录复制到纯英文路径 C:\Temp\mpnn_work\ 运行
然后把结果拷回项目目录
"""
import subprocess, json, sys, shutil
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
WORK = ROOT / "work" / "round4"
PROJECT_MPNN = WORK / "ProteinMPNN"

# === 用纯英文临时目录避开中文路径 bug ===
TMP_BASE = Path(r"C:\Temp\mpnn_work")
TMP_BASE.mkdir(parents=True, exist_ok=True)

TMP_MPNN = TMP_BASE / "ProteinMPNN"
TMP_PDB = TMP_BASE / "pdbs"
TMP_OUT = TMP_BASE / "out"

# 复制 (软链接更快, 但需管理员; 直接复制保险)
if not TMP_MPNN.exists():
    print("复制 ProteinMPNN 到 C:\\Temp\\mpnn_work...")
    shutil.copytree(PROJECT_MPNN, TMP_MPNN)
    print("  OK")

TMP_PDB.mkdir(exist_ok=True)
TMP_OUT.mkdir(exist_ok=True)
shutil.copy(WORK / "pdbs" / "2B3P.pdb", TMP_PDB / "2B3P.pdb")

helper = TMP_MPNN / "helper_scripts"
parsed = TMP_OUT / "parsed.jsonl"
chains = TMP_OUT / "chains.jsonl"
fixed = TMP_OUT / "fixed.jsonl"

# ============================================================
# Step 1: 解析 PDB
# ============================================================
print("Step 1: 解析 PDB...")
r = subprocess.run([
    sys.executable, str(helper / "parse_multiple_chains.py"),
    "--input_path", str(TMP_PDB),
    "--output_path", str(parsed),
], capture_output=True, text=True)
print(" return:", r.returncode)
if r.returncode != 0:
    print(" STDERR:", r.stderr); sys.exit(1)

# ============================================================
# Step 2: 链分配
# ============================================================
print("\nStep 2: 链分配...")
r = subprocess.run([
    sys.executable, str(helper / "assign_fixed_chains.py"),
    "--input_path", str(parsed),
    "--output_path", str(chains),
    "--chain_list", "A",
], capture_output=True, text=True)
print(" return:", r.returncode)

# ============================================================
# Step 3: 固定位置 (chromophore + 关键功能残基)
# ============================================================
# sfGFP 2B3P 实际起始残基: pos 1=M? PDB 序列可能从M开始
# chromophore TYG 在 pos 65-67 (1-based, 注意PDB编号可能与序列编号有偏移)
# 检查 parsed jsonl 看实际位置
with open(parsed) as f:
    p_data = json.loads(f.readline())
seq_a = p_data.get("seq_chain_A", "")
print(f"\n  2B3P chain A seq len = {len(seq_a)}")
# 找 TYG/SYG 位置
for cb in ["TYG", "SYG", "GYG"]:
    idx = seq_a.find(cb)
    if idx >= 0:
        print(f"  chromophore {cb} 在 1-based pos {idx+1}-{idx+3}")
        break

# 关键残基 (1-based, 基于序列)
# T65 Y66 G67 chromophore
# R96 chromophore 邻位
# H148 proton acceptor
# T203 chromophore 邻位
# E222 proton wire
chromo_pos = idx + 1 if idx >= 0 else 65
fixed_positions = " ".join(str(p) for p in [
    chromo_pos-4, chromo_pos-3, chromo_pos-2, chromo_pos-1,
    chromo_pos, chromo_pos+1, chromo_pos+2,
    chromo_pos+3, chromo_pos+4,
    96-4, 96, 148, 203, 222,
])
# 简化避免越界, 用绝对位置 (sfGFP)
fixed_positions = "61 62 63 64 65 66 67 68 69 70 96 145 148 167 203 205 222"
print(f"\n  固定位点 (1-based): {fixed_positions}")

print("\nStep 3: 生成 fixed_positions_dict...")
r = subprocess.run([
    sys.executable, str(helper / "make_fixed_positions_dict.py"),
    "--input_path", str(parsed),
    "--output_path", str(fixed),
    "--chain_list", "A",
    "--position_list", fixed_positions,
], capture_output=True, text=True)
print(" return:", r.returncode)

# ============================================================
# Step 4: 跑 ProteinMPNN
# ============================================================
for temp, n_seq, tag in [("0.1", 30, "T01"), ("0.3", 30, "T03")]:
    out_sub = TMP_OUT / tag
    out_sub.mkdir(exist_ok=True)
    print(f"\nStep 4: ProteinMPNN T={temp} ({n_seq} 序列)...")
    cmd = [
        sys.executable, str(TMP_MPNN / "protein_mpnn_run.py"),
        "--jsonl_path", str(parsed),
        "--chain_id_jsonl", str(chains),
        "--fixed_positions_jsonl", str(fixed),
        "--out_folder", str(out_sub),
        "--num_seq_per_target", str(n_seq),
        "--sampling_temp", temp,
        "--seed", "37",
        "--batch_size", "8",
        "--save_score", "1",
        "--model_name", "v_48_020",
        "--path_to_model_weights", str(TMP_MPNN / "vanilla_model_weights"),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    print(f"  return: {r.returncode}")
    if r.returncode != 0:
        print(" STDERR:", r.stderr[-500:])
        print(" STDOUT tail:", r.stdout[-500:])
    else:
        print(" OK")
        # 查找输出
        fa_files = list((out_sub / "seqs").glob("*.fa"))
        for fa in fa_files:
            with open(fa) as f:
                content = f.read()
            n_seqs = len([l for l in content.split("\n") if l.startswith(">")])
            print(f"  {fa.name}: {n_seqs} 个序列条目")

# ============================================================
# 拷贝结果回项目目录
# ============================================================
proj_out = WORK / "mpnn_output_v3"
if proj_out.exists():
    shutil.rmtree(proj_out)
shutil.copytree(TMP_OUT, proj_out)
print(f"\n✓ 结果已拷贝到 {proj_out}")
