"""
Round 8: LigandMPNN with chromophore (CUDA RTX 5080)
"""
import json, sys, time, os
from pathlib import Path
import torch

ROOT = Path(r"D:\生信\2026Protein Design")
R4 = ROOT / "work" / "round4"
R8 = ROOT / "work" / "round8"
R8.mkdir(parents=True, exist_ok=True)
LMPNN = ROOT / "work" / "round5" / "LigandMPNN"

print(f"CUDA: {torch.cuda.is_available()}")
sys.path.insert(0, str(LMPNN))
os.chdir(str(LMPNN))

# 固定残基: 只固定生色团关键位
FIXED = "A66 A96 A222"
CKPT = "./model_params/ligandmpnn_v_32_020_25.pt"

import shutil
shutil.copy(str(R4 / "pdbs" / "2B3P.pdb"), "2B3P.pdb")
print(f"Fixed: {FIXED} | Checkpoint: {os.path.getsize(CKPT)/1e6:.0f}MB")

from run import main as mpnn_main
import argparse

TASKS = [("T01", "0.1", 10), ("T03", "0.3", 10), ("T05", "0.5", 10)]

for tag, temp, n_batches in TASKS:
    tag_dir = Path(f"C:\\Temp\\r8_lmpnn_{tag}")
    tag_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n--- LigandMPNN T={temp} ---")
    t0 = time.time()
    
    args = argparse.Namespace(
        model_type="ligand_mpnn",
        checkpoint_ligand_mpnn=CKPT,
        pdb_path="2B3P.pdb",
        out_folder=str(tag_dir),
        temperature=float(temp),
        batch_size=5,
        number_of_batches=n_batches,
        fixed_residues=FIXED,
        omit_AA="C",
        seed=111,
        save_stats=0,
        verbose=0,
        ligand_mpnn_use_atom_context=1,
        ligand_mpnn_cutoff_for_score=8.0,
        ligand_mpnn_use_side_chain_context=0,
        pack_side_chains=0,
        num_cores=4,
        
        # Fill required defaults
        checkpoint_protein_mpnn="",
        checkpoint_per_residue_label_membrane_mpnn="",
        checkpoint_global_label_membrane_mpnn="",
        checkpoint_soluble_mpnn="",
        pdb_path_multi="",
        redesigned_residues="",
        chains_to_design="",
        parse_these_chains_only="",
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
        zero_indexed=0,
        file_ending="",
        checkpoint_path_sc="",
        number_of_packs_per_design=4,
        sc_num_denoising_steps=3,
        sc_num_samples=16,
        repack_everything=0,
        force_hetatm=0,
        packed_suffix="_packed",
        pack_with_ligand_context=1,
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
        print(f"  ❌ {str(e)[:200]}")

# Parse results
WT = "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

candidates, seen = [], set()
for tag, _, _ in TASKS:
    tag_dir = Path(f"C:\\Temp\\r8_lmpnn_{tag}")
    fa = tag_dir / "seqs" / "2B3P.fa"
    if not fa.exists():
        alt = list(tag_dir.rglob("*.fa"))
        fa = alt[0] if alt else None
    if not fa:
        continue
    with open(fa) as f:
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
            candidates.append({
                "name": f"R8_LMPNN_{tag}_{len(candidates)+1:03d}",
                "scaffold": "sfGFP_LMPNN",
                "n_muts": sum(1 for a,b in zip(full, WT) if a!=b),
                "length": len(full), "seq": full,
            })

print(f"\nCandidates: {len(candidates)}")
if candidates:
    candidates.sort(key=lambda x: x["n_muts"])
    print(f"Muts: {candidates[0]['n_muts']}-{candidates[-1]['n_muts']}")
    for c in candidates[:5]:
        print(f"  {c['name']}: {c['n_muts']} muts")

with open(R8 / "candidates_lmpnn_r8.json", "w", encoding="utf-8") as f:
    json.dump(candidates, f, indent=2, ensure_ascii=False)

os.chdir(str(ROOT))
print(f"Saved: {R8 / 'candidates_lmpnn_r8.json'}")
