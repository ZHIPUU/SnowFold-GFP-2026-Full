"""Round 5 P0-1: 加载扩展 LigandMPNN 候选 + 选 Top 准备 ESMFold"""
import re, json, pandas as pd
from pathlib import Path
import glob

ROOT = Path(r"D:\生信\2026Protein Design")
R5 = ROOT / "work" / "round5"

avGFP_wt = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
excl_seqs = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())

# 扫描 lmpnn_expanded
all_cands = []
exp_dir = R5 / "lmpnn_expanded"
for sub in sorted(exp_dir.glob("av_lmpnn_*")):
    fa_files = glob.glob(str(sub / "seqs" / "*.fa"))
    if not fa_files: continue
    with open(fa_files[0]) as f:
        content = f.read()
    blocks = [b for b in content.strip().split(">") if b.strip()]
    if not blocks: continue
    wt_seq = blocks[0].split("\n", 1)[1].strip().replace("\n", "")
    
    # offset 在 avGFP 中
    offset = avGFP_wt.find(wt_seq[:15])
    if offset < 0: offset = 2
    
    task_name = sub.name
    for i, b in enumerate(blocks[1:], 1):
        lines = b.split("\n", 1)
        if len(lines) < 2: continue
        header = lines[0].strip()
        seq = lines[1].strip().replace("\n", "")
        if len(seq) != len(wt_seq): continue
        if not any(cb in seq for cb in ["TYG","SYG","GYG"]): continue
        full_seq = avGFP_wt[:offset] + seq + avGFP_wt[offset + len(seq):]
        if len(full_seq) != 238: continue
        if not full_seq.startswith("M"): continue
        if set(full_seq) - set("ACDEFGHIKLMNPQRSTVWY"): continue
        if full_seq in excl_seqs: continue
        
        m_oc = re.search(r"overall_confidence=([\d.]+)", header)
        m_lc = re.search(r"ligand_confidence=([\d.]+)", header)
        m_rec = re.search(r"seq_rec=([\d.]+)", header)
        n_mut = sum(1 for a, b in zip(full_seq, avGFP_wt) if a != b)
        
        all_cands.append({
            "name": f"R5e_{task_name}_{i:03d}",
            "seq": full_seq,
            "scaffold": "avGFP_LMPNN",
            "n_muts": n_mut,
            "length": len(full_seq),
            "lmpnn_overall": float(m_oc.group(1)) if m_oc else None,
            "lmpnn_ligand": float(m_lc.group(1)) if m_lc else None,
            "lmpnn_recovery": float(m_rec.group(1)) if m_rec else None,
            "expected_tm": 80,
            "role": "lmpnn_de_novo",
            "task": task_name,
            "notes": f"LigandMPNN expanded {task_name}",
        })

print(f"总扩展候选: {len(all_cands)}")

# 去重 (seq)
seen = set()
dedup = []
for c in all_cands:
    if c["seq"] not in seen:
        seen.add(c["seq"])
        dedup.append(c)
print(f"去重后: {len(dedup)}")

# 按 ligand_confidence 取 top 30 进 ESMFold (避免太耗时)
dedup.sort(key=lambda x: -(x["lmpnn_ligand"] or 0))
top = dedup[:30]
print(f"\nTop 30 by lig_confidence:")
for c in top[:15]:
    print(f"  {c['name']:<35} mut={c['n_muts']:>3} lig={c['lmpnn_ligand']:.3f} oc={c['lmpnn_overall']:.3f}")

with open(R5 / "lmpnn_expanded_top30.json", "w", encoding="utf-8") as f:
    json.dump(top, f, indent=2, ensure_ascii=False)
print(f"\n保存 lmpnn_expanded_top30.json")
