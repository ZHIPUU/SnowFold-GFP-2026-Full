"""
Round 4 Step C-合并 v2: 修正 MPNN 输出对齐
==================================
关键发现:
  - 2B3P 实际 PDB chain A 序列从 pos 2 开始 (缺 M0)
  - MPNN 完全重设计了非固定位置 (包括 N 端)
  - 我们固定了 65-67 chromophore + R96 等关键位
  - MPNN 输出 = 设计后的全长 (231 aa, 缺 M)
  - 需要前缀 "M" 让其 ≥220 aa 且 M 起始
"""
import re, json, pandas as pd
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
WORK = ROOT / "work" / "round4"

sfgfp_wt = "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

# 读取 mpnn output 但保留 WT 序列做参考 (FA 文件第一条总是 WT)
def parse_fa(fa_path):
    """解析 ProteinMPNN .fa 输出, 返回 [(header, seq), ...]"""
    with open(fa_path) as f:
        content = f.read()
    entries = []
    blocks = content.strip().split(">")
    for block in blocks:
        if not block.strip(): continue
        lines = block.strip().split("\n", 1)
        if len(lines) < 2: continue
        header = lines[0].strip()
        seq = lines[1].strip().replace("\n", "")
        entries.append((header, seq))
    return entries


# 检查所有候选的 chromophore 位置
all_mpnn = []
for tag in ["T01", "T03", "T05"]:
    fa = WORK / "mpnn_output_final" / tag / "seqs" / "2B3P.fa"
    if not fa.exists(): continue
    entries = parse_fa(fa)
    
    # 第一条是 WT
    if entries:
        wt_header, wt_seq = entries[0]
        print(f"\n[{tag}] WT (PDB原序列): len={len(wt_seq)}, 前20={wt_seq[:20]}")
        # 找 chromophore
        for cb in ["TYG", "SYG", "GYG"]:
            if cb in wt_seq:
                pos = wt_seq.index(cb)
                print(f"  {cb} at 0-based pos {pos} (1-based {pos+1})")
                break
    
    # 后面的是设计序列
    for i, (header, seq) in enumerate(entries[1:], 1):
        m = re.search(r"score=([\d.]+)", header)
        score = float(m.group(1)) if m else None
        m2 = re.search(r"seq_recovery=([\d.]+)", header)
        recovery = float(m2.group(1)) if m2 else None
        
        # 验证 chromophore 保留
        has_chromo = any(c in seq for c in ["TYG", "SYG", "GYG"])
        
        all_mpnn.append({
            "tag": tag, "idx": i, "seq": seq, "len": len(seq),
            "score": score, "recovery": recovery, "has_chromo": has_chromo,
        })

print(f"\n总 MPNN 设计序列: {len(all_mpnn)}")
chromo_keep = [s for s in all_mpnn if s["has_chromo"]]
print(f"chromophore 保留的: {len(chromo_keep)}/{len(all_mpnn)}")

# ============================================================
# 加 M 前缀, 验证, 加入候选池
# ============================================================
excl_seqs = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())

candidates = []
for s in chromo_keep:
    # MPNN 输出已经 231 aa, 加 M = 232 aa (220-250 范围内)
    if s["seq"].startswith("M"):
        full_seq = s["seq"]
    else:
        full_seq = "M" + s["seq"]
    
    # 验证
    issues = []
    if not full_seq.startswith("M"): issues.append("no-M")
    if not (220 <= len(full_seq) <= 250): issues.append(f"len={len(full_seq)}")
    if set(full_seq) - set("ACDEFGHIKLMNPQRSTVWY"): issues.append("bad-AA")
    
    has_chromo = any(c in full_seq for c in ["TYG", "SYG", "GYG"])
    if not has_chromo: issues.append("no-chromo")
    
    in_excl = full_seq in excl_seqs
    
    # 与 sfGFP 的差异 (无意义但记录)
    # MPNN 完全重设计了非固定位, 与 sfGFP WT 差异巨大 (~60 mut)
    
    if not issues and not in_excl:
        candidates.append({
            "name": f"MPNN_{s['tag']}_{s['idx']:03d}",
            "seq": full_seq,
            "scaffold": "sfGFP_MPNN",
            "role": "de_novo_MPNN",
            "n_muts": -1,  # 无意义, de novo 设计
            "length": len(full_seq),
            "expected_tm": 80,
            "mpnn_score": s["score"],
            "mpnn_recovery": s["recovery"],
            "mpnn_tag": s["tag"],
            "notes": f"ProteinMPNN T={s['tag']} recovery={s['recovery']:.3f} de novo",
        })

print(f"\n通过验证: {len(candidates)} / {len(chromo_keep)}")
print(f"  长度分布: {set(c['length'] for c in candidates)}")

# 按 recovery 降序排, 取 top 15 进入 ESMFold 评估
candidates.sort(key=lambda x: -(x["mpnn_recovery"] or 0))
top_candidates = candidates[:15]  # 评估太多耗时, 选 top 15
print(f"\nTop 15 (按 recovery):")
for c in top_candidates:
    print(f"  {c['name']:<22} recovery={c['mpnn_recovery']:.3f} score={c['mpnn_score']:.3f}")

# 保存
with open(WORK / "mpnn_candidates.json", "w", encoding="utf-8") as f:
    json.dump(top_candidates, f, indent=2, ensure_ascii=False)
print(f"\n✓ 保存 {len(top_candidates)} 条到 mpnn_candidates.json")
