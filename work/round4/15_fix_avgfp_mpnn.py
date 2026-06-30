"""
Round 4: 检查 PDB 真实序列位置 + 重跑 avGFP 用正确的固定位
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
# 检查 avGFP 2WUR 解析后的序列
# ============================================================
av_dir = TMP_BASE / "avGFP_2WUR"
with open(av_dir / "parsed.jsonl") as f:
    p = json.loads(f.readline())
seq = p["seq_chain_A"]
print(f"2WUR seq len = {len(seq)}")
print(f"前 50: {seq[:50]}")
for cb in ["TYG", "SYG", "GYG"]:
    if cb in seq:
        idx = seq.index(cb)
        print(f"{cb} at 1-based pos {idx+1} (0-based {idx})")

# 2WUR 序列 chromophore 实际位置
# 修正固定位
chromo_idx = next((seq.index(cb) for cb in ["TYG","SYG","GYG"] if cb in seq), -1)
print(f"\nchromophore 0-based pos = {chromo_idx}")
chromo_1based = chromo_idx + 1
# 固定 chromophore ± 5 + 几个关键位
fixed_positions = []
for offset in range(-5, 6):
    fixed_positions.append(chromo_1based + offset)
# 添加关键功能位 (相对 chromophore 偏移)
# R96 ~ chromophore+31, E222 ~ chromophore+157
fixed_positions.extend([
    chromo_1based + 31,   # R96 equivalent
    chromo_1based + 80,   # H148 equivalent
    chromo_1based + 138,  # T203 equivalent
    chromo_1based + 157,  # E222 equivalent
])
fixed_positions = sorted(set([p for p in fixed_positions if 1 <= p <= len(seq)]))
fixed_str = " ".join(str(p) for p in fixed_positions)
print(f"修正后固定位: {fixed_str}")

# 重新生成 fixed.jsonl
r = subprocess.run([
    sys.executable, str(helper / "make_fixed_positions_dict.py"),
    "--input_path", "parsed.jsonl", "--output_path", "fixed_v2.jsonl",
    "--chain_list", "A", "--position_list", fixed_str,
], capture_output=True, text=True, cwd=str(av_dir))
print(f"重新生成 fixed_v2: return={r.returncode}")

# 重新跑
for temp, n_seq, tag in [("0.1", 20, "T01_v2"), ("0.3", 15, "T03_v2")]:
    out_sub = f"out_{tag}"
    (av_dir / out_sub).mkdir(exist_ok=True)
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
    ], capture_output=True, text=True, timeout=600, cwd=str(av_dir))
    if r.returncode != 0:
        print(" FAIL:", r.stderr[-300:])
    else:
        fa_files = list((av_dir / out_sub / "seqs").glob("*.fa"))
        for fa in fa_files:
            n = sum(1 for l in open(fa).read().split("\n") if l.startswith(">"))
            print(f"  ✓ {fa.name}: {n} 序列")

# 拷贝到项目目录
proj_out = WORK / "mpnn_multi_scaffold" / "avGFP_2WUR_v2"
if proj_out.exists():
    shutil.rmtree(proj_out)
proj_out.mkdir()
for tag in ["T01_v2", "T03_v2"]:
    src = av_dir / f"out_{tag}"
    if src.exists():
        shutil.copytree(src, proj_out / tag)

# ============================================================
# 同样检查 mBaoJin (8QBJ)
# ============================================================
print("\n" + "="*60)
print("检查 mBaoJin (8QBJ)")
print("="*60)
mb_dir = TMP_BASE / "mBaoJin_8QBJ"
with open(mb_dir / "parsed.jsonl") as f:
    p = json.loads(f.readline())
seq = p["seq_chain_A"]
print(f"8QBJ seq len = {len(seq)}")
print(f"前 30: {seq[:30]}")
for cb in ["TYG", "SYG", "GYG"]:
    if cb in seq:
        idx = seq.index(cb)
        print(f"{cb} at 1-based pos {idx+1}")
        chromo_idx = idx

print(f"\n现有固定位: 63 64 65 66 67 68 69 70 96 140 141 165 171 201")
print(f"chromophore 0-based = {chromo_idx}")
# 默认 mBaoJin 序列 GYG 在 pos 66-68, 那 1-based 66-68
# 之前固定的 63-68 包含了 chromophore
print(f"  原固定位 63-68 是否覆盖 chromophore {chromo_idx+1}-{chromo_idx+3}? ", end="")
print("YES" if 63 <= chromo_idx+1 <= 68 else "NO")

print("\n✓ 完成")
