"""
Round 7: Fixed-residue ProteinMPNN Design + ESMFold Round-trip
==============================================================
核心策略 (基于文献):
  1. Sumida et al. JACS 2024: 固定关键残基 → pLDDT 从 50 → 85+
  2. EnhancedMPNN (ICLR 2026): 用 pLDDT 做 DPO 奖励信号, 3× 成功率
  
改进 vs Round 4:
  - 使用 `--use_soluble_model` (GFP 是可溶蛋白)
  - 固定更多残基: 生色团 + 口袋 + β-barrel 核心
  - 更多温度: T=0.1, 0.3, 0.5, 0.7
  - 更多序列: 50/温度
  - 纯英文路径 (C:\Temp\mpnn_r7)

 Pipeline:
   1. ProteinMPNN 生成 200 条新序列
   2. ESMFold round-trip 验证
   3. 严格过滤: pTM>0.75, pLDDT>80, 生色团>85
   4. 亮度估计 + Top-6 选择
"""
import subprocess, json, sys, shutil, re, time, numpy as np
from pathlib import Path
import pandas as pd

ROOT = Path(r"D:\生信\2026Protein Design")
R4 = ROOT / "work" / "round4"
OUT = ROOT / "work" / "round7"

# ============================================================
# 1. 准备 ProteinMPNN (英文路径)
# ============================================================
TMP_BASE = Path(r"C:\Temp\mpnn_r7")
TMP_MPNN = TMP_BASE / "ProteinMPNN"

if not TMP_MPNN.exists():
    print("复制 ProteinMPNN 到英文路径...")
    shutil.copytree(R4 / "ProteinMPNN", TMP_MPNN)
    print("  OK")

helper = TMP_MPNN / "helper_scripts"
PDB_SRC = R4 / "pdbs" / "2B3P.pdb"
assert PDB_SRC.exists(), f"PDB 不存在: {PDB_SRC}"

# ============================================================
# 2. 解析 PDB + 准备固定残基
# ============================================================
# sfGFP (2B3P) 残基编号 (1-based)
# 生色团核心区域 (规则要求 58-72, 210-230):
# 固定策略:
#   Tier 1 (必须固定): 生色团 + 成熟催化 + 质子传递
#     Y66, G67 (chromophore 三联体), R96 (maturation), E222 (proton wire)
#   Tier 2 (强烈推荐固定): 生色团邻近环境
#     T65, L64, Q69, S72, H148, T203 (chromophore 结合袋)
#   Tier 3 (β-barrel 核心): sfGFP 折叠核心
#     S30R, F64L, F99S, M153T, V163A, I171V, A206V

OUT.mkdir(parents=True, exist_ok=True)

# 解析 PDB 链
parsed_jsonl = TMP_BASE / "parsed.jsonl"
print("\nStep 1: 解析 PDB 链...")
r = subprocess.run([
    sys.executable, str(helper / "parse_multiple_chains.py"),
    "--input_path", str(PDB_SRC.parent),
    "--output_path", str(parsed_jsonl),
], capture_output=True, text=True)
print(f"  returncode: {r.returncode}")
if r.stderr:
    print(f"  stderr: {r.stderr[:300]}")

# 修复: 替换中文路径名为简单名称 (ProteinMPNN 输出文件名会用到)
with open(parsed_jsonl) as f:
    parsed_data = json.loads(f.readline())
# 改名为简单英文名
parsed_data["name"] = "2B3P"
with open(parsed_jsonl, "w") as f:
    json.dump(parsed_data, f)
seq_len = len(parsed_data["seq_chain_A"])
print(f"  PDB 序列长度: {seq_len}")

# ============================================================
# 3. 准备固定位置文件
# ============================================================
# sfGFP 关键残基 (1-based 编号)
# 固定 Tier 1 + 2 + 3:
# 生色团区域 (58-72): 尽量多固定
# 生色团结对袋 (210-230): 固定关键位
# R96, E222 等
FIXED_POSITIONS = [
    # Tier 1: 生色团 + 催化
    65, 66, 67,   # TYG chromophore
    96,           # R96 chromophore maturation
    222,          # E222 proton wire
    # Tier 2: 生色团环境
    64,           # F/L64 靠近 chromophore
    69,           # Q69 chromophore 环境
    72,           # S/A72 chromophore 环境
    148,          # H148 质子受体
    145,          # Y/F145 生色团邻位
    167,          # T167 氢键网络
    203,          # T203 β-barrel
    205,          # S205 核心
    # Tier 3: β-barrel 核心折叠
    30,           # S/R30 关键稳定突变
    99,           # F/S99 折叠核心
    153,          # M/T153 折叠核心
    163,          # V/A163 折叠核心
    171,          # I/V171 稳定
    206,          # A/V206 防二聚
]

# 过滤超范围的位置
FIXED_POSITIONS = sorted([p for p in FIXED_POSITIONS if 1 <= p <= seq_len])
fixed_str = " ".join(str(p) for p in FIXED_POSITIONS)
print(f"\nStep 2: 固定 {len(FIXED_POSITIONS)} 个残基:")
print(f"  Positions: {fixed_str}")

# 生成 fixed_positions_dict.jsonl
# 注意: chain_id_jsonl 需要为空(单链场景)或用正确格式
# ProteinMPNN 空字符串=不加载 chain_id, 默认所有链都可见面
fixed_jsonl = TMP_BASE / "fixed.jsonl"

r = subprocess.run([
    sys.executable, str(helper / "make_fixed_positions_dict.py"),
    "--input_path", str(parsed_jsonl),
    "--output_path", str(fixed_jsonl),
    "--chain_list", "A",
    "--position_list", fixed_str,
], capture_output=True, text=True)
print(f"  fixed dict: returncode={r.returncode}")

# ============================================================
# 4. 运行 ProteinMPNN (可溶性模型)
# ============================================================
# 温度: 0.1(保守), 0.3, 0.5(平衡), 0.7(激进)
TEMPS = [
    ("0.1", 50, "T01"),
    ("0.3", 50, "T03"),
    ("0.5", 50, "T05"),
    ("0.7", 50, "T07"),
]

all_fa_files = []

for temp, n_seq, tag in TEMPS:
    out_dir = TMP_BASE / tag
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\nStep 3: ProteinMPNN T={temp} ({n_seq} seqs, 可溶性模型)...")
    t0 = time.time()
    
    # 显式指定权重路径 (Windows 路径修复)
    model_weights = str(TMP_MPNN / "soluble_model_weights")
    
    cmd = [
        sys.executable, str(TMP_MPNN / "protein_mpnn_run.py"),
        "--jsonl_path", str(parsed_jsonl),
        "--chain_id_jsonl", "",
        "--fixed_positions_jsonl", str(fixed_jsonl),
        "--out_folder", str(out_dir),
        "--num_seq_per_target", str(n_seq),
        "--sampling_temp", temp,
        "--seed", "42",
        "--batch_size", "10",
        "--save_score", "1",
        "--model_name", "v_48_020",
        "--path_to_model_weights", model_weights,
    ]
    
    print(f"  cmd: {' '.join(str(c) for c in cmd[:6])} ...")
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=1200)
    elapsed = time.time() - t0
    
    if r.returncode != 0:
        print(f"  ❌ 失败 (rc={r.returncode}), stderr: {r.stderr[-500:]}")
    else:
        print(f"  ✅ 完成 ({elapsed:.0f}s)")
    
    # 收集输出的 .fa 文件
    fa_path = out_dir / "seqs" / "2B3P.fa"
    if fa_path.exists():
        all_fa_files.append(fa_path)
        print(f"  输出: {fa_path}")

print(f"\n{'='*80}")
print(f"ProteinMPNN 设计完成! 共 {len(all_fa_files)} 个输出文件")
print(f"{'='*80}")

# ============================================================
# 5. 解析 ProteinMPNN 输出序列 + 填回固定位 AA
# ============================================================
# ProteinMPNN 在固定位置输出 X, 需要用 WT 序列填回
# 2B3P PDB 实际从 sfGFP 的 S2 开始 (缺少 N 端 M)
# PDB 序列 (231 aa): SKGEELF... 对应 sfGFP 的 pos 2-238
# 因此 PDB pos 1 = sfGFP pos 2 (S)
# 完整 sfGFP WT (238 aa): MSKGEELF...
SFGFP_WT = "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
# PDB 起始偏移 (PDB pos 1 = WT pos 1, 但 PDB 序列是 231 aa 从 S 开始)
# 实际 PDB 序列对应 WT pos 2-238 (跳过 N 端 M)
PDB_OFFSET = 1  # PDB seq[0] = WT seq[1] = 'S'

def fill_x_and_align(seq_pdb, wt_seq):
    """
    将 MPNN 输出的 PDB-长度序列 (231 aa, 从 S 开始) 对齐到完整 WT (238 aa, 从 M 开始)
    并填回 X 占位符.
    
    MPNN 输出: 231 aa (S...), 固定位为 X
    目标: 238 aa (MSKGEEL...), 固定位用 WT 填回
    """
    # 第一步: 填回 X (用 PDB 对应位置的 WT AA)
    seq_filled = list(seq_pdb)
    for i, aa in enumerate(seq_filled):
        if aa == "X":
            # PDB pos i 对应 WT pos i+1 (因为 PDB 从 S2 开始)
            wt_idx = i + PDB_OFFSET
            if wt_idx < len(wt_seq):
                seq_filled[i] = wt_seq[wt_idx]
    
    # 第二步: 补 N 端 M, 补回完整长度
    full_seq = "M" + "".join(seq_filled)
    
    # 如果长度不够, 从 WT 补足
    if len(full_seq) < 220:
        full_seq = full_seq + wt_seq[len(full_seq):220]
    
    return full_seq[:250]  # 截断到 250

def count_muts(seq, wt_seq):
    """计数相对于 WT 的突变数 (对齐比较)"""
    min_len = min(len(seq), len(wt_seq))
    return sum(1 for i in range(min_len) if seq[i] != wt_seq[i])

candidates = []
seen_seqs = set()

for fa_path in all_fa_files:
    tag = fa_path.parent.parent.name  # T01, T03, etc.
    with open(fa_path) as f:
        content = f.read()
    
    blocks = content.strip().split(">")
    for i, block in enumerate(blocks):
        if not block.strip():
            continue
        lines = block.strip().split("\n", 1)
        if len(lines) < 2:
            continue
        header = lines[0].strip()
        seq = lines[1].strip().replace("\n", "")
        
        # 跳过 WT 序列
        if "WT" in header.upper() or i == 0:
            continue
        
        # 填回 X + 对齐到完整 WT 长度
        seq = fill_x_and_align(seq, SFGFP_WT)
        
        # 合规检查
        if seq[0] != "M" or len(seq) < 220 or len(seq) > 250:
            continue
        if set(seq) - set("ACDEFGHIKLMNPQRSTVWY"):
            continue
        
        if seq not in seen_seqs:
            seen_seqs.add(seq)
            n_muts = count_muts(seq, SFGFP_WT)
            
            # 简单过滤: 突变数太少的可能是 WT 本身, 太多则风险高
            if n_muts < 5:
                continue  # 太接近 WT
            if n_muts > 200:
                continue  # 可能对齐错误
            candidates.append({
                "name": f"R7_MPNN_{tag}_{i:03d}",
                "scaffold": "sfGFP_MPNN",
                "n_muts": n_muts,
                "length": len(seq),
                "seq": seq,
                "notes": f"ProteinMPNN soluble T={tag[1:]}.{tag[2]} fixed={len(FIXED_POSITIONS)}res",
            })

print(f"\n解析完成: {len(candidates)} 条新候选 (去重)")

# ============================================================
# 6. ESMFold 验证
# ============================================================
# 注意: esm v3/v2 冲突, 改用 subprocess 调用 esmfold 的 transformers API
# 实际上 ESMFold 需要用到 fair-esm 的 API, 而当前环境中 esm v3 冲突
# 回退方案: 用我们已经存在的候选评估方式

print(f"\n{'='*80}")
print(f"Round 7: 保存候选, ESMFold 需后续用 transformers API 单独跑")
print(f"共 {len(candidates)} 条候选")
print(f"{'='*80}")

# 保存候选
cand_path = OUT / "candidates_round7.json"
with open(cand_path, "w", encoding="utf-8") as f:
    json.dump(candidates, f, indent=2, ensure_ascii=False)
print(f"候选已保存: {cand_path}")

# 统计
n_muts_list = [c["n_muts"] for c in candidates]
print(f"\n突变数统计:")
print(f"  平均: {np.mean(n_muts_list):.0f}")
print(f"  范围: {min(n_muts_list)} - {max(n_muts_list)}")
print(f"  <50: {sum(1 for m in n_muts_list if m < 50)}")
print(f"  <80: {sum(1 for m in n_muts_list if m < 80)}")

# 找出最接近 WT 的 5 条 (突变最少)
candidates.sort(key=lambda x: x["n_muts"])
print(f"\n突变最少 Top-5:")
for c in candidates[:5]:
    print(f"  {c['name']}: {c['n_muts']} muts, len={c['length']}")

print(f"\n⚠️ 重要提示: ESMFold 验证需要修复 esm v3/v2 冲突")
print(f"  临时方案: pip install fair-esm (已安装)")
print(f"  或手动运行: transformers API 加载 facebook/esmfold_v1")
print(f"\n下一步: 修复 ESMFold 后运行 ESMFold round-trip 验证")
