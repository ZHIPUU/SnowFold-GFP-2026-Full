"""
Round 4: 重跑 mBaoJin MPNN 用正确 chromophore 位置 + 提取所有 MPNN 结果
"""
import subprocess, json, sys, shutil
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
WORK = ROOT / "work" / "round4"

TMP_BASE = Path(r"C:\Temp\mpnn_multi")
TMP_MPNN = TMP_BASE / "ProteinMPNN"
helper = TMP_MPNN / "helper_scripts"
weights_path = str(TMP_MPNN / "vanilla_model_weights")

# ============================================================
# 重跑 mBaoJin 用正确的 chromophore pos 36-38
# ============================================================
mb_dir = TMP_BASE / "mBaoJin_8QBJ"
with open(mb_dir / "parsed.jsonl") as f:
    p = json.loads(f.readline())
seq = p["seq_chain_A"]

chromo_idx = next(seq.index(cb) for cb in ["TYG","SYG","GYG"] if cb in seq)
chromo_1 = chromo_idx + 1
print(f"mBaoJin GYG 在 1-based pos {chromo_1}-{chromo_1+2}")

# 固定 chromophore ± 5 + 关键功能位 (相对偏移)
fixed_pos = []
for off in range(-5, 6):
    fixed_pos.append(chromo_1 + off)
# StayGold 家族的关键: R134 (chromophore +98 类似), E222 等价
fixed_pos.extend([chromo_1 + 60, chromo_1 + 100, chromo_1 + 130, chromo_1 + 160])
fixed_pos = sorted(set([p for p in fixed_pos if 1 <= p <= len(seq)]))
fixed_str = " ".join(str(p) for p in fixed_pos)
print(f"修正后固定位: {fixed_str}")

r = subprocess.run([
    sys.executable, str(helper / "make_fixed_positions_dict.py"),
    "--input_path", "parsed.jsonl", "--output_path", "fixed_v2.jsonl",
    "--chain_list", "A", "--position_list", fixed_str,
], capture_output=True, text=True, cwd=str(mb_dir))
print(f"fixed_v2: return={r.returncode}")

for temp, n_seq, tag in [("0.1", 20, "T01_v2"), ("0.3", 15, "T03_v2")]:
    out_sub = f"out_{tag}"
    (mb_dir / out_sub).mkdir(exist_ok=True)
    print(f"\nT={temp} ({n_seq} 序列)...")
    r = subprocess.run([
        sys.executable, str(TMP_MPNN / "protein_mpnn_run.py"),
        "--jsonl_path", "parsed.jsonl",
        "--chain_id_jsonl", "chains.jsonl",
        "--fixed_positions_jsonl", "fixed_v2.jsonl",
        "--out_folder", out_sub,
        "--num_seq_per_target", str(n_seq),
        "--sampling_temp", temp,
        "--seed", "37",
        "--batch_size", "8",
        "--save_score", "1",
        "--model_name", "v_48_020",
        "--path_to_model_weights", weights_path,
    ], capture_output=True, text=True, timeout=600, cwd=str(mb_dir))
    if r.returncode == 0:
        fa = mb_dir / out_sub / "seqs" / "input.fa"
        if fa.exists():
            n = sum(1 for l in open(fa).read().split("\n") if l.startswith(">"))
            print(f"  ✓ {n} 序列")
    else:
        print(" FAIL:", r.stderr[-200:])

# 拷贝到项目目录
proj_out = WORK / "mpnn_multi_scaffold" / "mBaoJin_8QBJ_v2"
if proj_out.exists():
    shutil.rmtree(proj_out)
proj_out.mkdir(parents=True)
for tag in ["T01_v2", "T03_v2"]:
    src = mb_dir / f"out_{tag}"
    if src.exists():
        shutil.copytree(src, proj_out / tag)
print(f"\n✓ mBaoJin v2 拷贝到 {proj_out}")
