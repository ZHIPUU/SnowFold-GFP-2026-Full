"""
Round 8: ProteinMPNN direct API (CUDA + RTX 5080)
直接加载 checkpoint，绕开 run.py 的中文路径编码问题
"""
import json, sys, time, os, numpy as np
from pathlib import Path
import torch

ROOT = Path(r"D:\生信\2026Protein Design")
R4 = ROOT / "work" / "round4"
R8 = ROOT / "work" / "round8"
R8.mkdir(parents=True, exist_ok=True)
R5 = ROOT / "work" / "round5"

print(f"CUDA: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"GPU: {torch.cuda.get_device_name(0)}")

# 切换到 LigandMPNN 目录 (为了导入)
sys.path.insert(0, str(R5 / "LigandMPNN"))
os.chdir(str(R5 / "LigandMPNN"))

from run import main as mpnn_main
import argparse

# 固定残基: 生色团 + 催化 + 核心骨架 (17个)
FIXED = "A58 A65 A66 A67 A68 A69 A70 A71 A72 A96 A148 A171 A203 A205 A206 A219 A222"
print(f"固定残基: {FIXED}")

# 复制 PDB
import shutil
shutil.copy(str(R4 / "pdbs" / "2B3P.pdb"), "2B3P.pdb")

TASKS = [("T01", "0.1", 10), ("T02", "0.2", 10), ("T03", "0.3", 10)]

for tag, temp, n_batches in TASKS:
    tag_dir = Path(f"C:\\Temp\\r8_mpnn_{tag}")
    tag_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n--- ProteinMPNN T={temp} ---")
    t0 = time.time()
    
    # 构建 args
    args = argparse.Namespace(
        model_type="protein_mpnn",
        checkpoint_protein_mpnn="./model_params/proteinmpnn_v_48_020.pt",
        checkpoint_ligand_mpnn="./model_params/ligandmpnn_v_32_010_25.pt",
        checkpoint_per_residue_label_membrane_mpnn="./model_params/per_residue_label_membrane_mpnn_v_48_020.pt",
        checkpoint_global_label_membrane_mpnn="./model_params/global_label_membrane_mpnn_v_48_020.pt",
        checkpoint_soluble_mpnn="./model_params/solublempnn_v_48_020.pt",
        pdb_path="2B3P.pdb",
        pdb_path_multi="",
        fixed_residues=FIXED,
        redesigned_residues="",
        omit_AA="C",
        chains_to_design="",
        parse_these_chains_only="",
        seed=111,
        temperature=float(temp),
        batch_size=5,
        number_of_batches=n_batches,
        out_folder=str(tag_dir),
        save_stats=0,
        zero_indexed=0,
        file_ending="",
        verbose=1,
        fasta_seq_separation=":",
        bias_AA="",
        bias_AA_per_residue="",
        homo_oligomer=0,
        symmetry_residues="",
        symmetry_weights="",
        transmembrane_buried="",
        transmembrane_interface="",
        global_transmembrane_label=0,
        parse_atoms_with_zero_occupancy=0,
        pack_side_chains=0,
        checkpoint_path_sc="./model_params/ligandmpnn_sc_v_32_002_16.pt",
        number_of_packs_per_design=4,
        sc_num_denoising_steps=3,
        sc_num_samples=16,
        repack_everything=0,
        force_hetatm=0,
        packed_suffix="_packed",
        pack_with_ligand_context=1,
        ligand_mpnn_use_atom_context=1,
        ligand_mpnn_cutoff_for_score=8.0,
        ligand_mpnn_use_side_chain_context=0,
        fixed_residues_multi="",
        redesigned_residues_multi="",
        bias_AA_per_residue_multi="",
        omit_AA_per_residue="",
        omit_AA_per_residue_multi="",
    )
    
    try:
        mpnn_main(args)
        print(f"  ✅ {time.time()-t0:.0f}s")
    except Exception as e:
        print(f"  ❌ {e}")

# 解析序列
WT = "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

candidates = []
seen = set()

for tag, _, _ in TASKS:
    tag_dir = Path(f"C:\\Temp\\r8_mpnn_{tag}")
    fa_path = tag_dir / "seqs" / "2B3P.fa"
    if not fa_path.exists():
        alt = list(tag_dir.rglob("*.fa"))
        fa_path = alt[0] if alt else None
    if not fa_path:
        continue
    
    with open(fa_path) as f:
        content = f.read()
    
    blocks = content.strip().split(">")
    for block in blocks:
        if not block.strip():
            continue
        lines = block.strip().split("\n", 1)
        if len(lines) < 2:
            continue
        seq = lines[1].strip().replace(" ", "").replace("\n", "")
        
        seq_list = list(seq)
        for j, aa in enumerate(seq_list):
            if aa == "X" and j + 1 < len(WT):
                seq_list[j] = WT[j + 1]
            elif aa == "X":
                seq_list[j] = "A"
        full_seq = "M" + "".join(seq_list)
        
        if full_seq[0] != "M" or len(full_seq) < 220 or len(full_seq) > 250:
            continue
        if full_seq not in seen:
            seen.add(full_seq)
            n_muts = sum(1 for a, b in zip(full_seq, WT) if a != b)
            candidates.append({
                "name": f"R8_PMPNN_{tag}_{len(candidates)+1:03d}",
                "scaffold": "sfGFP_MPNN",
                "n_muts": n_muts,
                "length": len(full_seq),
                "seq": full_seq,
                "notes": f"ProteinMPNN direct T={tag.replace('T','')}"
            })

print(f"\n解析: {len(candidates)} 条")
if candidates:
    candidates.sort(key=lambda x: x["n_muts"])
    print(f"突变: {candidates[0]['n_muts']}-{candidates[-1]['n_muts']}")
    for c in candidates[:5]:
        print(f"  {c['name']}: {c['n_muts']} muts")

with open(R8 / "candidates_mpnn_r8.json", "w", encoding="utf-8") as f:
    json.dump(candidates, f, indent=2, ensure_ascii=False)

os.chdir(str(ROOT))
