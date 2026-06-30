"""
Round 4 Step C-合并 v3: 用 PDB WT 填回 X 占位符
==================================
关键发现:
  - 2B3P 实际序列从 'S' 开始 (1-based pos 1=S), chromophore TYG 在 pos 64-66
  - 我们固定的 1-based 65-67 实际 = WT pos 65-67 = "YGV" (偏移1!)
  - ProteinMPNN 输出 X 占位符表示我们指定的"固定位置"
  - 需要用 PDB WT 序列在对应位置填回真实 AA
"""
import re, json, pandas as pd
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
WORK = ROOT / "work" / "round4"

sfgfp_wt = "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"


def parse_fa(fa_path):
    with open(fa_path) as f:
        content = f.read()
    entries = []
    for block in content.strip().split(">"):
        if not block.strip(): continue
        lines = block.strip().split("\n", 1)
        if len(lines) < 2: continue
        entries.append((lines[0].strip(), lines[1].strip().replace("\n", "")))
    return entries


def fill_x_with_wt(mpnn_seq, pdb_wt):
    """用 PDB WT 序列填回 MPNN 输出中的 X"""
    assert len(mpnn_seq) == len(pdb_wt), f"长度不匹配 {len(mpnn_seq)} vs {len(pdb_wt)}"
    return "".join(pdb_wt[i] if c == "X" else c for i, c in enumerate(mpnn_seq))


excl_seqs = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())

all_candidates = []
for tag in ["T01", "T03", "T05"]:
    fa = WORK / "mpnn_output_final" / tag / "seqs" / "2B3P.fa"
    if not fa.exists(): continue
    entries = parse_fa(fa)
    if not entries: continue
    
    # 第一条是 PDB WT (含 X 占位符)
    pdb_wt_with_x = entries[0][1]
    print(f"\n[{tag}] PDB WT (含X): {pdb_wt_with_x[:75]}...")
    print(f"  X 数量: {pdb_wt_with_x.count('X')}, 总长 {len(pdb_wt_with_x)}")
    
    # 但第一条 WT 也有 X (固定位置), 需要用真实 sfGFP 填回
    # sfGFP_wt 比 pdb_wt 多 N 端 M, 所以 sfGFP_wt[1:] 与 pdb_wt 对齐
    actual_wt = sfgfp_wt[1:]  # 去掉 M, 长度 237
    
    # pdb_wt 长度 231, sfGFP[1:] 长度 237, 差 6 → pdb 缺 N 端 5 个 + C 端 1 个
    # 实际看下: pdb_wt[:20] vs actual_wt[:20]
    print(f"  pdb_wt[:20]: {pdb_wt_with_x[:20]}")
    print(f"  sfGFP[1:21]: {actual_wt[:20]}")
    print(f"  sfGFP[5:25]: {actual_wt[5:25]}")
    # 实际对齐
    # 用模式匹配找对齐
    matches = []
    for offset in range(0, 15):
        # 比对前 30 个非 X 字符
        cmp_a = pdb_wt_with_x[:30]
        cmp_b = actual_wt[offset:offset+30]
        m = sum(1 for a, b in zip(cmp_a, cmp_b) if a == b or a == 'X')
        matches.append((offset, m))
    best_offset = max(matches, key=lambda x: x[1])
    print(f"  最佳对齐 offset = {best_offset[0]}, 匹配 {best_offset[1]}/30")
    
    # sfGFP[offset+1:] 对齐 pdb_wt
    align_start = best_offset[0] + 1  # +1 因为 actual_wt 已经去了 M
    wt_aligned = sfgfp_wt[align_start:align_start + len(pdb_wt_with_x)]
    
    # 后续设计
    for i, (header, mpnn_seq) in enumerate(entries[1:], 1):
        if len(mpnn_seq) != len(pdb_wt_with_x):
            continue
        
        # 用 sfGFP WT 填回 X
        filled = "".join(wt_aligned[j] if c == "X" else c for j, c in enumerate(mpnn_seq))
        
        # 加 M 前缀, 加 sfGFP 缺失的 N 端
        # 实际策略: 直接用 sfGFP 的前 align_start 字符 + filled
        full_seq = sfgfp_wt[:align_start] + filled + sfgfp_wt[align_start + len(filled):]
        # 现在 full_seq 应该是 238 aa, 与 sfGFP 等长
        
        # 验证
        if len(full_seq) != 238: continue
        if not full_seq.startswith("M"): continue
        if "X" in full_seq: continue
        if any(c not in "ACDEFGHIKLMNPQRSTVWY" for c in full_seq): continue
        if not any(cb in full_seq for cb in ["TYG", "SYG", "GYG"]): continue
        if full_seq in excl_seqs: continue
        
        m = re.search(r"score=([\d.]+)", header)
        score = float(m.group(1)) if m else None
        m2 = re.search(r"seq_recovery=([\d.]+)", header)
        recovery = float(m2.group(1)) if m2 else None
        
        # 突变数 (vs sfGFP WT)
        n_mut = sum(1 for a, b in zip(full_seq, sfgfp_wt) if a != b)
        
        all_candidates.append({
            "name": f"MPNN_{tag}_{i:03d}",
            "seq": full_seq,
            "scaffold": "sfGFP_MPNN",
            "role": "de_novo_MPNN",
            "n_muts": n_mut,
            "length": len(full_seq),
            "expected_tm": 80,
            "mpnn_score": score,
            "mpnn_recovery": recovery,
            "mpnn_tag": tag,
            "notes": f"ProteinMPNN T={tag} recovery={recovery:.3f} score={score:.3f}",
        })

print(f"\n\n总有效候选: {len(all_candidates)}")
if all_candidates:
    muts = [c["n_muts"] for c in all_candidates]
    print(f"突变数: min={min(muts)} max={max(muts)} mean={sum(muts)/len(muts):.1f}")

# 选 Top 15 (按 recovery 降序, 突变数 ≤ 100 防止 30%阈值风险)
filtered = [c for c in all_candidates if c["n_muts"] <= 80]
filtered.sort(key=lambda x: -(x["mpnn_recovery"] or 0))
top = filtered[:12]

print(f"\nTop 12 (突变数≤80, 按 recovery 降序):")
for c in top:
    print(f"  {c['name']:<22} n_mut={c['n_muts']:>3} recovery={c['mpnn_recovery']:.3f} score={c['mpnn_score']:.3f}")

with open(WORK / "mpnn_candidates_final.json", "w", encoding="utf-8") as f:
    json.dump(top, f, indent=2, ensure_ascii=False)
print(f"\n✓ 保存 {len(top)} 条到 mpnn_candidates_final.json")
