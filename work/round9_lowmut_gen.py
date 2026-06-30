"""
Round 9: Generate low-mutation (<30) variants of MPNN_T01_014
=============================================================
策略: 从 WT PDB 出发，用 --redesigned_residues 只设计 25 个位点
其余位点保持 WT → 总突变数 ≤25

挑选 25 个位点: 从 MPNN_T01_014 的 57 突变中选非生色团、表面位点
"""
import json, sys, time, os, random, shutil
from pathlib import Path
import torch
import numpy as np

ROOT = Path(r"D:\生信\2026Protein Design")
R4 = ROOT / "work" / "round4"
R9 = ROOT / "work" / "round9"
R9.mkdir(parents=True, exist_ok=True)
LMPNN = ROOT / "work" / "round5" / "LigandMPNN"

WT = "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

# 1. 获取 MPNN_T01_014 序列
with open(R9.parent / "round6" / "final_6_new_rules.json", encoding="utf-8") as f:
    data = json.load(f)
seq1 = [d["seq"] for d in data if "MPNN_T01_014" in d["name"]][0]

# 2. 找出 57 个突变位点 (PDB numbering: skip M1, so PDB pos = seq index)
mut_positions = []
for i, (a, b) in enumerate(zip(WT, seq1)):
    if a != b:
        mut_positions.append(i + 1)  # 0-based seq → 1-based PDB (M=1)
print(f"Total mutations in MPNN_T01_014: {len(mut_positions)}")

# 3. 排除生色团区域 (58-72) 和 核心 (210-230)
chromo_region = set(range(58, 73)) | set(range(210, 231))
filtered = [p for p in mut_positions if p not in chromo_region]
print(f"After chromo filter: {len(filtered)}")

# 4. 固定种子，保证可重复
random.seed(42)
pick_count = 25
# 选择均匀分布的位点 (按位置排序后等间距采样)
filtered.sort()
step = max(1, len(filtered) // pick_count)
selected = [filtered[i] for i in range(0, len(filtered), step)][:pick_count]
selected.sort()
print(f"Selected {len(selected)} positions for redesign:")
print(f"  {selected}")

# 5. 使用 --fixed_residues: 固定除选中的25个位点外的所有位点
all_positions = set(range(2, 239))  # PDB 2-238 (total 237)
fixed_positions = sorted(all_positions - set(selected))
fixed_str = " ".join(f"A{p}" for p in fixed_positions)
print(f"\nFixed {len(fixed_positions)} positions, designing {len(selected)} positions")
print(f"Fixed: {fixed_str[:100]}...")

# 6. 运行 LigandMPNN
print(f"\n{'='*60}")
print(f"CUDA: {torch.cuda.is_available()}")
sys.path.insert(0, str(LMPNN))
os.chdir(str(LMPNN))
shutil.copy(str(R4 / "pdbs" / "2B3P.pdb"), "2B3P.pdb")

CKPT = "./model_params/solublempnn_v_48_020.pt"
CKPT_TYPE = "soluble_mpnn"
TEMP, BATCH, N_BATCH = "0.2", 5, 10  # 50 seqs

tag_dir = Path("C:\\Temp\\r9_lowmut")
tag_dir.mkdir(parents=True, exist_ok=True)

from run import main as mpnn_main
import argparse

args = argparse.Namespace(
    model_type=CKPT_TYPE,
    checkpoint_soluble_mpnn=CKPT,
    pdb_path="2B3P.pdb",
    out_folder=str(tag_dir),
    temperature=float(TEMP),
    batch_size=BATCH,
    number_of_batches=N_BATCH,
    redesigned_residues="",
    fixed_residues=fixed_str,  # fix all except 25 design positions
    omit_AA="C",
    seed=42,
    save_stats=0, verbose=0,
    checkpoint_protein_mpnn="", checkpoint_per_residue_label_membrane_mpnn="",
    checkpoint_global_label_membrane_mpnn="",
    pdb_path_multi="",
    chains_to_design="",
    parse_these_chains_only="", fasta_seq_separation=":",
    bias_AA="", bias_AA_per_residue="", homo_oligomer=0,
    symmetry_residues="", symmetry_weights="",
    transmembrane_buried="", transmembrane_interface="",
    global_transmembrane_label=0, parse_atoms_with_zero_occupancy=0,
    zero_indexed=0, file_ending="",
    checkpoint_path_sc="", number_of_packs_per_design=4,
    sc_num_denoising_steps=3, sc_num_samples=16,
    repack_everything=0, force_hetatm=0, packed_suffix="_packed",
    pack_with_ligand_context=1,
    fixed_residues_multi="", redesigned_residues_multi="",
    bias_AA_per_residue_multi="", omit_AA_per_residue="",
    omit_AA_per_residue_multi="",
    num_cores=4,
)

t0 = time.time()
mpnn_main(args)
print(f"Design done: {time.time()-t0:.0f}s")

# 7. 解析序列
fa_path = tag_dir / "seqs" / "2B3P.fa"
if not fa_path.exists():
    alt = list(tag_dir.rglob("*.fa"))
    fa_path = alt[0] if alt else None

candidates = []
seen = set()
if fa_path:
    with open(fa_path) as f:
        content = f.read()
    for block in content.strip().split(">")[1:]:
        lines = block.strip().split("\n", 1)
        if len(lines) < 2: continue
        seq = lines[1].strip().replace(" ", "").replace("\n", "")
        seq_l = list(seq)
        for j, aa in enumerate(seq_l):
            if aa == "X" and j+1 < len(WT):
                seq_l[j] = WT[j+1]
        full = "M" + "".join(seq_l)
        if full[0] != "M" or len(full) < 220 or len(full) > 250:
            continue
        if full not in seen:
            seen.add(full)
            nmuts = sum(1 for a, b in zip(full, WT) if a != b)
            candidates.append({
                "name": f"R9_lowmut_{len(candidates)+1:03d}",
                "scaffold": "sfGFP_lowmut",
                "n_muts": nmuts, "length": len(full), "seq": full,
                "notes": f"fixed: {fixed_str}"
            })

# 过滤: 只保留突变数 < 30
candidates = [c for c in candidates if c["n_muts"] < 30]
print(f"\nParsed: {len(candidates)} candidates (<30 muts)")
for c in candidates[:5]:
    print(f"  {c['name']}: {c['n_muts']} muts")

# 8. 保存
with open(R9 / "candidates_r9.json", "w", encoding="utf-8") as f:
    json.dump(candidates, f, indent=2, ensure_ascii=False)
print(f"\nSaved: {R9 / 'candidates_r9.json'}")

os.chdir(str(ROOT))
