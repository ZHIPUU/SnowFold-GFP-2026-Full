"""
Round 4 优化继续: 跑 mBaoJin (8QBJ) + avGFP (2WUR) 的 ProteinMPNN 设计
==================================
策略:
  1. mBaoJin (8QBJ): 用真实 PDB 做 inverse folding, 可能 pLDDT 从 39 -> 60+
  2. avGFP (2WUR): 多骨架多样性
  3. 固定 chromophore + 关键功能位点
"""
import subprocess, json, sys, shutil, os
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
WORK = ROOT / "work" / "round4"

TMP_BASE = Path(r"C:\Temp\mpnn_multi")
TMP_BASE.mkdir(parents=True, exist_ok=True)

TMP_MPNN = TMP_BASE / "ProteinMPNN"
if not TMP_MPNN.exists():
    print("复制 ProteinMPNN...")
    shutil.copytree(WORK / "ProteinMPNN", TMP_MPNN)
    print("  OK")

helper = TMP_MPNN / "helper_scripts"
weights_path = str(TMP_MPNN / "vanilla_model_weights")

# ============================================================
# 任务定义
# ============================================================
tasks = [
    {
        "name": "mBaoJin_8QBJ",
        "pdb_src": WORK / "pdbs" / "8QBJ.pdb",
        # mBaoJin chromophore GYG 在 pos 66-68 (基于序列扫描)
        # 关键: R96 chromophore 邻位, E222 类似
        # 固定 chromophore + 周围 5 个 + Q140P, H141Q (单体化关键), C165Y, N171Y, T201A
        # 注意 mBaoJin 编号与 sfGFP 不同
        "fixed": "63 64 65 66 67 68 69 70 96 140 141 165 171 201",
    },
    {
        "name": "avGFP_2WUR",
        "pdb_src": WORK / "pdbs" / "2WUR.pdb",
        # avGFP chromophore TYG/SYG 在 pos 65-67
        "fixed": "61 62 63 64 65 66 67 68 69 70 96 145 148 167 203 205 222",
    },
]

all_designs = []

for task in tasks:
    print(f"\n{'='*60}")
    print(f"任务: {task['name']}")
    print('='*60)
    
    task_dir = TMP_BASE / task["name"]
    task_dir.mkdir(exist_ok=True)
    
    # 复制 PDB
    pdb_dst = task_dir / "input.pdb"
    shutil.copy(task["pdb_src"], pdb_dst)
    
    # Step 1: 解析
    print("Step 1: 解析 PDB...")
    r = subprocess.run([
        sys.executable, str(helper / "parse_multiple_chains.py"),
        "--input_path", ".", "--output_path", "parsed.jsonl",
    ], capture_output=True, text=True, cwd=str(task_dir))
    if r.returncode != 0:
        print(" FAIL:", r.stderr[-300:]); continue
    
    with open(task_dir / "parsed.jsonl") as f:
        p = json.loads(f.readline())
    seq_a = p.get("seq_chain_A", "")
    print(f"  序列长度 = {len(seq_a)}")
    for cb in ["TYG", "SYG", "GYG"]:
        if cb in seq_a:
            print(f"  {cb} at 1-based pos {seq_a.index(cb)+1}")
            break
    
    # Step 2: 链分配
    r = subprocess.run([
        sys.executable, str(helper / "assign_fixed_chains.py"),
        "--input_path", "parsed.jsonl", "--output_path", "chains.jsonl",
        "--chain_list", "A",
    ], capture_output=True, text=True, cwd=str(task_dir))
    
    # Step 3: 固定位置
    print(f"Step 3: 固定位点 {task['fixed']}")
    r = subprocess.run([
        sys.executable, str(helper / "make_fixed_positions_dict.py"),
        "--input_path", "parsed.jsonl", "--output_path", "fixed.jsonl",
        "--chain_list", "A", "--position_list", task["fixed"],
    ], capture_output=True, text=True, cwd=str(task_dir))
    
    # Step 4: 跑 ProteinMPNN
    for temp, n_seq, tag in [("0.1", 20, "T01"), ("0.3", 15, "T03")]:
        out_sub = f"out_{tag}"
        (task_dir / out_sub).mkdir(exist_ok=True)
        print(f"\nStep 4: T={temp} ({n_seq} 序列)...")
        cmd = [
            sys.executable, str(TMP_MPNN / "protein_mpnn_run.py"),
            "--jsonl_path", "parsed.jsonl",
            "--chain_id_jsonl", "chains.jsonl",
            "--fixed_positions_jsonl", "fixed.jsonl",
            "--out_folder", out_sub,
            "--num_seq_per_target", str(n_seq),
            "--sampling_temp", temp,
            "--seed", "37",
            "--batch_size", "8",
            "--save_score", "1",
            "--model_name", "v_48_020",
            "--path_to_model_weights", weights_path,
        ]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=600, cwd=str(task_dir))
        if r.returncode != 0:
            print(" FAIL:", r.stderr[-300:])
        else:
            fa_files = list((task_dir / out_sub / "seqs").glob("*.fa"))
            for fa in fa_files:
                with open(fa) as f:
                    content = f.read()
                n_seqs = sum(1 for l in content.split("\n") if l.startswith(">"))
                print(f"  ✓ {fa.name}: {n_seqs} 序列")

# ============================================================
# 拷贝结果回项目目录
# ============================================================
proj_out = WORK / "mpnn_multi_scaffold"
if proj_out.exists():
    shutil.rmtree(proj_out)
proj_out.mkdir()
for task in tasks:
    src = TMP_BASE / task["name"]
    if src.exists():
        # 只拷贝输出
        dst_task = proj_out / task["name"]
        dst_task.mkdir(exist_ok=True)
        for tag in ["T01", "T03"]:
            src_out = src / f"out_{tag}"
            if src_out.exists():
                shutil.copytree(src_out, dst_task / tag)
print(f"\n✓ 结果拷贝到 {proj_out}")
