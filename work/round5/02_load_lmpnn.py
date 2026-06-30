"""
Round 5 P0-1b: 解析 LigandMPNN 输出 + 排除列表过滤
注意:
  - LigandMPNN 输出的是完整序列 (含 chromophore)
  - 长度比 sfGFP WT (238) 略短 (223 = 缺前7 + 后8), 需要补
  - 前缀 sfGFP[0:7] = "MSKGEEL", 后缀 sfGFP[-8:] = "GMDELYK" 等
"""
import re, json, pandas as pd
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
R5 = ROOT / "work" / "round5"

# WT 序列
sfGFP_wt = "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
avGFP_wt = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

print(f"sfGFP WT length: {len(sfGFP_wt)}")
print(f"avGFP WT length: {len(avGFP_wt)}")

# ============================================================
# 解析 LigandMPNN .fa
# ============================================================
def parse_fa(path):
    """解析 .fa, 返回 [(seq, conf, lig_conf, recovery), ...]"""
    with open(path) as f:
        content = f.read()
    entries = []
    blocks = content.strip().split(">")
    for block in blocks:
        if not block.strip(): continue
        lines = block.strip().split("\n", 1)
        if len(lines) < 2: continue
        header = lines[0].strip()
        seq = lines[1].strip().replace("\n", "")
        # 提取置信度
        m_oc = re.search(r"overall_confidence=([\d.]+)", header)
        m_lc = re.search(r"ligand_confidence=([\d.]+)", header)
        m_rec = re.search(r"seq_rec=([\d.]+)", header)
        oc = float(m_oc.group(1)) if m_oc else None
        lc = float(m_lc.group(1)) if m_lc else None
        rec = float(m_rec.group(1)) if m_rec else None
        is_wt = oc is None  # 第一条是 WT, 无置信度
        entries.append({"seq": seq, "oc": oc, "lc": lc, "rec": rec, "is_wt": is_wt})
    return entries


# 加载所有 LigandMPNN 输出
tasks = {
    "sfGFP_lmpnn_T01": ("sfGFP", "2B3P.fa"),
    "sfGFP_lmpnn_T03": ("sfGFP", "2B3P.fa"),
    "avGFP_lmpnn_T01": ("avGFP", "2WUR.fa"),
}

excl_seqs = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())

all_cands = []
for task_name, (scaffold, fa_name) in tasks.items():
    fa = R5 / "lmpnn_output" / task_name / "seqs" / fa_name
    if not fa.exists():
        print(f"  ✗ {task_name}: 文件不存在")
        continue
    entries = parse_fa(fa)
    print(f"\n[{task_name}] {len(entries)} 条 (含 WT)")
    if not entries: continue
    
    # 第一条是 WT (含 chromophore 但缺 N/C 端的可能)
    wt_entry = entries[0]
    wt_seq = wt_entry["seq"]
    print(f"  PDB WT seq length: {len(wt_seq)} (前 30: {wt_seq[:30]})")
    
    # 找 PDB WT 在 WT 完整序列中的 offset
    ref_wt = sfGFP_wt if scaffold == "sfGFP" else avGFP_wt
    
    # 用前 20 个字符 (假设 LigandMPNN 不改 N 端区) 找 offset
    offset = ref_wt.find(wt_seq[:20])
    if offset < 0:
        # 第一条可能也被设计过, 找不到完全匹配, 用 PDB 中 chromophore SER65 锚定
        # 暂时假设 offset = 2 (类似 Round 4)
        offset = 2
    print(f"  在 {scaffold} 中 offset = {offset}")
    
    # 处理设计序列
    for i, e in enumerate(entries[1:], 1):
        if e["is_wt"]: continue
        if len(e["seq"]) != len(wt_seq):
            continue
        # 补 N/C 端
        full_seq = ref_wt[:offset] + e["seq"] + ref_wt[offset + len(e["seq"]):]
        # 检查
        if len(full_seq) != len(ref_wt): continue
        if not full_seq.startswith("M"): continue
        if not (220 <= len(full_seq) <= 250): continue
        if set(full_seq) - set("ACDEFGHIKLMNPQRSTVWY"): continue
        if not any(cb in full_seq for cb in ["TYG", "SYG", "GYG"]): continue
        if full_seq in excl_seqs: continue
        
        n_mut = sum(1 for a, b in zip(full_seq, ref_wt) if a != b)
        all_cands.append({
            "name": f"R5_{task_name}_{i:03d}",
            "seq": full_seq,
            "scaffold": scaffold + "_LMPNN",
            "n_muts": n_mut,
            "length": len(full_seq),
            "lmpnn_overall_confidence": e["oc"],
            "lmpnn_ligand_confidence": e["lc"],
            "lmpnn_recovery": e["rec"],
            "task": task_name,
            "expected_tm": 80,  # 后续 ThermoMPNN 调整
            "role": "lmpnn_de_novo",
            "notes": f"LigandMPNN {task_name} oc={e['oc']:.3f} lc={e['lc']:.3f} rec={e['rec']:.3f}",
        })

print(f"\n\n总通过验证: {len(all_cands)}")

# 按 ligand_confidence 排序 (LigandMPNN 独有的 chromophore-aware 评分)
all_cands.sort(key=lambda x: -(x["lmpnn_ligand_confidence"] or 0))

print(f"\nTop 15 (按 ligand_confidence):")
for c in all_cands[:15]:
    print(f"  {c['name']:<35} n_mut={c['n_muts']:>3} lig_conf={c['lmpnn_ligand_confidence']:.3f} oc={c['lmpnn_overall_confidence']:.3f}")

# 保留 top 20 进入 ESMFold 验证
top = all_cands[:20]
with open(R5 / "lmpnn_candidates.json", "w", encoding="utf-8") as f:
    json.dump(top, f, indent=2, ensure_ascii=False)

print(f"\n保存 top {len(top)} 到 lmpnn_candidates.json")
