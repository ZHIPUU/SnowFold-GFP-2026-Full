"""
Round 5 P0-1c: 加载 avGFP LigandMPNN v2 (chromophore 100% 保留) + ESMFold 评估
"""
import re, json, pandas as pd
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
R5 = ROOT / "work" / "round5"

# avGFP WT
avGFP_wt = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

# 加载 LigandMPNN v2 avGFP 结果
fa_path = R5 / "lmpnn_output_v2" / "avGFP_lmpnn_v2_T01" / "seqs" / "2WUR.fa"
with open(fa_path) as f:
    content = f.read()
blocks = [b for b in content.strip().split(">") if b.strip()]
wt_seq = blocks[0].split("\n", 1)[1].strip().replace("\n", "")
print(f"2WUR PDB seq len = {len(wt_seq)}")
print(f"2WUR 前 20: {wt_seq[:20]}")

# 找 PDB WT 在 avGFP WT 中的 offset
offset = avGFP_wt.find(wt_seq[:15])
if offset < 0:
    # 用前 5 个字符
    offset = avGFP_wt.find(wt_seq[:5])
print(f"avGFP offset = {offset}")

excl_seqs = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())

cands = []
for i, b in enumerate(blocks[1:], 1):
    lines = b.split("\n", 1)
    if len(lines) < 2: continue
    header = lines[0].strip()
    seq = lines[1].strip().replace("\n", "")
    
    if len(seq) != len(wt_seq): continue
    if not any(cb in seq for cb in ["TYG","SYG","GYG"]): continue
    
    # 补 N/C 端
    full_seq = avGFP_wt[:offset] + seq + avGFP_wt[offset + len(seq):]
    if len(full_seq) != 238: continue
    if not full_seq.startswith("M"): continue
    if set(full_seq) - set("ACDEFGHIKLMNPQRSTVWY"): continue
    if full_seq in excl_seqs: continue
    
    # 提取置信度
    m_oc = re.search(r"overall_confidence=([\d.]+)", header)
    m_lc = re.search(r"ligand_confidence=([\d.]+)", header)
    m_rec = re.search(r"seq_rec=([\d.]+)", header)
    
    n_mut = sum(1 for a, b in zip(full_seq, avGFP_wt) if a != b)
    
    cands.append({
        "name": f"R5_av_lmpnn_v2_{i:03d}",
        "seq": full_seq,
        "scaffold": "avGFP_LMPNN",
        "n_muts": n_mut,
        "length": len(full_seq),
        "lmpnn_overall": float(m_oc.group(1)) if m_oc else None,
        "lmpnn_ligand": float(m_lc.group(1)) if m_lc else None,
        "lmpnn_recovery": float(m_rec.group(1)) if m_rec else None,
        "expected_tm": 80,
        "role": "lmpnn_de_novo",
        "notes": f"LigandMPNN v2 avGFP chromo-aware",
    })

print(f"\n通过验证: {len(cands)}")
if cands:
    cands.sort(key=lambda x: -(x["lmpnn_ligand"] or 0))
    print(f"\nTop 10 by ligand_confidence:")
    for c in cands[:10]:
        print(f"  {c['name']:<28} n_mut={c['n_muts']:>3} lig={c['lmpnn_ligand']:.3f} oc={c['lmpnn_overall']:.3f} rec={c['lmpnn_recovery']:.3f}")

# 取 top 15 进入 ESMFold
top = cands[:15]
with open(R5 / "lmpnn_v2_candidates.json", "w", encoding="utf-8") as f:
    json.dump(top, f, indent=2, ensure_ascii=False)
print(f"\n保存 {len(top)} 条到 lmpnn_v2_candidates.json")
