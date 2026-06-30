"""
Round 4 优化: 解析 avGFP MPNN 输出 + ESMFold 评估
==================================
注意:
  - 2WUR chain A 起始 pos = 36 (chromophore TYG)
  - 2WUR 起始残基对应 avGFP WT pos 1+offset, 需要确定 offset
"""
import json, re, pandas as pd
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
WORK = ROOT / "work" / "round4"

# 加载 avGFP WT
avGFP_wt = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
print(f"avGFP WT len: {len(avGFP_wt)}")

# 2WUR 解析后的实际序列
av_seq = "KGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGK..."  # 从前 50 知道
# 在 avGFP_wt 中找对齐
print(f"avGFP[2:52]: {avGFP_wt[2:52]}")
# 看起来 2WUR 缺 前 2 个 MS, 从 K (pos 3) 开始

# 解析 2WUR MPNN 输出
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

all_av_mpnn = []
for tag in ["T01_v2", "T03_v2"]:
    fa = WORK / "mpnn_multi_scaffold" / "avGFP_2WUR_v2" / tag / "seqs" / "input.fa"
    if not fa.exists():
        print(f"  {tag} 不存在")
        continue
    entries = parse_fa(fa)
    if not entries: continue
    
    # 第一条是 WT (含 X 占位符)
    wt_header, wt_with_x = entries[0]
    print(f"\n[{tag}] WT len={len(wt_with_x)}, X count={wt_with_x.count('X')}")
    print(f"  前 50: {wt_with_x[:50]}")
    
    # 验证对齐: 2WUR 序列与 avGFP_wt[2:] 对齐
    # 找 chromophore TYG/SYG 在 wt_with_x 中位置 (应该是 X)
    # 取 WT 的非 X 部分, 与 avGFP 比对
    # 实际 av_seq[:35] 应该匹配 avGFP_wt[2:37]
    
    # 用 avGFP_wt 填回 X
    # offset: avGFP_wt 中 2WUR seq[0]='K' 对应的位置
    avgfp_offset = avGFP_wt.find(wt_with_x[:10].replace("X", ""))
    if avgfp_offset < 0:
        # 找前 5 个非 X 字符
        non_x = "".join(c for c in wt_with_x[:15] if c != "X")[:10]
        avgfp_offset = avGFP_wt.find(non_x)
    print(f"  在 avGFP WT 中 offset = {avgfp_offset}")
    
    # avGFP WT[offset:] 对齐 2WUR seq
    wt_aligned = avGFP_wt[avgfp_offset:avgfp_offset + len(wt_with_x)]
    
    # 验证对齐
    matches = sum(1 for a, b in zip(wt_with_x, wt_aligned) if a == b or a == 'X')
    print(f"  对齐匹配: {matches}/{len(wt_with_x)}")
    
    # 处理设计序列
    for i, (header, mpnn_seq) in enumerate(entries[1:], 1):
        if len(mpnn_seq) != len(wt_with_x):
            continue
        # 用 avGFP WT 填回 X
        filled = "".join(wt_aligned[j] if c == "X" else c for j, c in enumerate(mpnn_seq))
        # 补回 N 端缺失
        full_seq = avGFP_wt[:avgfp_offset] + filled + avGFP_wt[avgfp_offset + len(filled):]
        
        # 验证
        if len(full_seq) != 238: continue
        if not full_seq.startswith("M"): continue
        if "X" in full_seq: continue
        if not any(cb in full_seq for cb in ["TYG","SYG","GYG"]): continue
        if set(full_seq) - set("ACDEFGHIKLMNPQRSTVWY"): continue
        
        m = re.search(r"score=([\d.]+)", header)
        score = float(m.group(1)) if m else None
        m2 = re.search(r"seq_recovery=([\d.]+)", header)
        recovery = float(m2.group(1)) if m2 else None
        
        n_mut = sum(1 for a, b in zip(full_seq, avGFP_wt) if a != b)
        
        all_av_mpnn.append({
            "name": f"MPNN_av_{tag}_{i:03d}",
            "seq": full_seq,
            "scaffold": "avGFP_MPNN",
            "role": "de_novo_MPNN",
            "n_muts": n_mut,
            "length": len(full_seq),
            "expected_tm": 80,
            "mpnn_score": score,
            "mpnn_recovery": recovery,
            "mpnn_tag": tag,
            "notes": f"ProteinMPNN avGFP T={tag} rec={recovery:.3f}",
        })

print(f"\n总 avGFP MPNN 有效候选: {len(all_av_mpnn)}")

# 检查不在排除列表
excl_seqs = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())
all_av_mpnn = [c for c in all_av_mpnn if c["seq"] not in excl_seqs]
print(f"通过排除列表: {len(all_av_mpnn)}")

# 按 recovery 排序, 取 top 8
all_av_mpnn.sort(key=lambda x: -(x["mpnn_recovery"] or 0))
top_av = [c for c in all_av_mpnn if c["n_muts"] <= 80][:8]

print(f"\nTop 8 avGFP MPNN:")
for c in top_av:
    print(f"  {c['name']:<28} n_mut={c['n_muts']:>3} rec={c['mpnn_recovery']:.3f}")

with open(WORK / "mpnn_avgfp_candidates.json", "w", encoding="utf-8") as f:
    json.dump(top_av, f, indent=2, ensure_ascii=False)
print(f"\n保存到 mpnn_avgfp_candidates.json")
