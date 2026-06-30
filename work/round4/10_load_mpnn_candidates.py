"""
Round 4 Step C-合并: 解析 ProteinMPNN 输出, 加入候选池
==================================
注意:
  - ProteinMPNN 输出的序列长度可能与 2B3P 晶体结构相同 (231 aa, 缺 N 端 7个)
  - 需要补上 N 端 "MSKGEEL" 等才能 ≥220 aa
  - 提交要求 220-250 aa, M 开头
"""
import re, json, pandas as pd
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
WORK = ROOT / "work" / "round4"

# ============================================================
# 解析 ProteinMPNN 输出的 .fa 文件
# ============================================================
sfgfp_wt = "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
print(f"sfGFP WT 长度: {len(sfgfp_wt)}")

mpnn_seqs = []
for tag in ["T01", "T03", "T05"]:
    fa = WORK / "mpnn_output_final" / tag / "seqs" / "2B3P.fa"
    if not fa.exists():
        continue
    with open(fa) as f:
        content = f.read()
    
    blocks = content.strip().split(">")
    for i, block in enumerate(blocks):
        if not block.strip(): continue
        lines = block.strip().split("\n", 1)
        if len(lines) < 2: continue
        header = lines[0].strip()
        seq = lines[1].strip().replace("\n", "")
        # 跳过 header 第一条是 WT
        if i == 1 or "score=" not in header:
            continue
        # 提取 score
        m = re.search(r"score=([\d.]+)", header)
        score = float(m.group(1)) if m else None
        m2 = re.search(r"seq_recovery=([\d.]+)", header)
        recovery = float(m2.group(1)) if m2 else None
        
        mpnn_seqs.append({
            "tag": tag,
            "idx": i,
            "seq_short": seq,  # 231 aa, 缺 N 端
            "score": score,
            "recovery": recovery,
            "len_short": len(seq),
        })

print(f"\n总共加载 ProteinMPNN 序列: {len(mpnn_seqs)}")
if mpnn_seqs:
    print("Top 5 (按 recovery 排序):")
    for s in sorted(mpnn_seqs, key=lambda x: -(x['recovery'] or 0))[:5]:
        print(f"  [{s['tag']}] recovery={s['recovery']:.3f} score={s['score']:.3f} len={s['len_short']}")

# ============================================================
# 补 N 端 → 238 aa (sfGFP长度)
# ============================================================
# sfGFP WT N 端前 7 aa = "MSKGEEL"
# 2B3P 晶体可能从 pos 8 开始, 或是 pos 1 但 N 端缺失
# 验证: 比对 mpnn 序列 vs sfGFP WT, 找对齐位置
def find_align_offset(mpnn_seq, wt):
    """找 mpnn_seq 在 wt 中的最佳对齐起点"""
    best_offset = 0
    best_match = 0
    for offset in range(0, len(wt) - len(mpnn_seq) + 1):
        match = sum(1 for i in range(len(mpnn_seq)) if mpnn_seq[i] == wt[offset+i])
        if match > best_match:
            best_match = match
            best_offset = offset
    return best_offset, best_match

if mpnn_seqs:
    sample = mpnn_seqs[0]["seq_short"]
    offset, n_match = find_align_offset(sample, sfgfp_wt)
    print(f"\nSample MPNN seq 在 sfGFP WT 中对齐: offset={offset}, match={n_match}/{len(sample)}")

# 假设 2B3P 起始 pos 对应 sfGFP WT pos 8 (N端缺7) 或 1
# 实际看下样本: mpnn_seq[0:5] vs sfGFP[0:5]
print(f"  MPNN seq 起始: '{mpnn_seqs[0]['seq_short'][:10]}'")
print(f"  sfGFP[0:10]: '{sfgfp_wt[:10]}'")
print(f"  sfGFP[7:17]: '{sfgfp_wt[7:17]}'")

# 决定: 默认 MPNN 序列对应 sfGFP[7:] (即 pos 8-238)
# 在前面补回 "MSKGEEL" (7 aa)
n_term = sfgfp_wt[:7]  # MSKGEEL

# ============================================================
# 构建候选 + 验证
# ============================================================
excl_seqs = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())

candidates = []
for i, s in enumerate(mpnn_seqs):
    # 补 N 端
    full_seq = n_term + s["seq_short"]  # 7 + 231 = 238 aa
    
    # 验证
    issues = []
    if not full_seq.startswith("M"): issues.append("no-M")
    if not (220 <= len(full_seq) <= 250): issues.append(f"len={len(full_seq)}")
    if set(full_seq) - set("ACDEFGHIKLMNPQRSTVWY"): issues.append("bad-AA")
    
    # chromophore
    has_chromo = any(c in full_seq for c in ["TYG", "SYG", "GYG"])
    if not has_chromo: issues.append("no-chromo")
    
    in_excl = full_seq in excl_seqs
    
    # 与 sfGFP WT 突变数
    n_mut = sum(1 for a, b in zip(full_seq, sfgfp_wt) if a != b)
    
    if not issues and not in_excl:
        candidates.append({
            "name": f"MPNN_{s['tag']}_{s['idx']:03d}",
            "seq": full_seq,
            "scaffold": "sfGFP_MPNN",
            "role": "de_novo_MPNN",
            "n_muts": n_mut,
            "length": len(full_seq),
            "expected_tm": 80,  # ProteinMPNN 通常保守, 设为 sfGFP baseline
            "mpnn_score": s["score"],
            "mpnn_recovery": s["recovery"],
            "mpnn_tag": s["tag"],
            "notes": f"ProteinMPNN T={s['tag']} recovery={s['recovery']:.3f} score={s['score']:.3f}",
        })

print(f"\n通过验证: {len(candidates)} / {len(mpnn_seqs)} (有 chromophore + 不在排除列表)")

if candidates:
    print("\n突变数分布:")
    from collections import Counter
    muts = [c["n_muts"] for c in candidates]
    print(f"  min={min(muts)} max={max(muts)} mean={sum(muts)/len(muts):.1f}")
    
    # 按突变数升序排, recovery 高的优先
    candidates.sort(key=lambda x: (x["n_muts"], -x["mpnn_recovery"]))
    print("\nTop 10 (按突变数升序):")
    for c in candidates[:10]:
        print(f"  {c['name']:<22} n_mut={c['n_muts']:>3} recovery={c['mpnn_recovery']:.3f}")
    
    # 保存
    with open(WORK / "mpnn_candidates.json", "w", encoding="utf-8") as f:
        json.dump(candidates, f, indent=2, ensure_ascii=False)
    print(f"\n保存 {len(candidates)} 条到 mpnn_candidates.json")
