"""全量排除列表检查 + 提交合规验证"""
import pandas as pd
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")

# 加载数据
excl = pd.read_csv(ROOT / "Exclusion_List.csv")
excl_seqs = set(excl["Sequence"].astype(str).str.strip())
print(f"排除列表: {len(excl)} 条, 唯一: {len(excl_seqs)} 条")

sub = pd.read_csv(ROOT / "submission_yourteamname.csv")
print(f"提交候选: {len(sub)} 条\n")

# 检查
print("=" * 60)
print("排除列表检查")
print("=" * 60)
for _, row in sub.iterrows():
    seq = str(row["Sequence"]).strip()
    in_excl = seq in excl_seqs
    status = "❌ 命中排除列表!" if in_excl else "✓ 通过"
    print(f"Seq {row['Seq_ID']}: {status} (长度={len(seq)})")

print("\n" + "=" * 60)
print("合规性验证 (长度220-250, M开头, 标准AA, chromophore)")
print("=" * 60)
for _, row in sub.iterrows():
    seq = str(row["Sequence"]).strip()
    name = row.get("Team_Name", row.get("Seq_ID", ""))
    issues = []
    if not seq.startswith("M"):
        issues.append("不以M开头")
    if len(seq) < 220 or len(seq) > 250:
        issues.append(f"长度{len(seq)}不在220-250")
    invalid_aa = {c for c in seq if c not in "ACDEFGHIKLMNPQRSTVWY"}
    if invalid_aa:
        issues.append(f"含非标准AA: {invalid_aa}")
    # chromophore check
    chromo = None
    for tri in ["TYG", "SYG", "GYG", "CYG", "HYG"]:
        if tri in seq:
            chromo = tri
            break
    if not chromo:
        issues.append("chromophore三联体缺失")

    if issues:
        print(f"Seq {row['Seq_ID']}: {' | '.join(issues)}")
    else:
        print(f"Seq {row['Seq_ID']}: ✓ 全部合规 (chromophore={chromo})")

# 统计每个 scaffold 的突变数
print("\n" + "=" * 60)
print("当前提交概览")
print("=" * 60)

# 加载WT序列
with open(ROOT / "AAseqs of 5 GFP proteins_20260511.txt") as f:
    wt_text = f.read()

wts = {}
for block in wt_text.split(">"):
    block = block.strip()
    if not block:
        continue
    lines = block.split("\n")
    name = lines[0].strip()
    seq = "".join(l.strip() for l in lines[1:] if l.strip() and not l.strip().startswith("#"))
    wts[name] = seq

for _, row in sub.iterrows():
    seq = str(row["Sequence"]).strip()
    # 找最匹配的WT
    best_wt, best_match = None, 0
    for wt_name, wt_seq in wts.items():
        match = sum(1 for a, b in zip(seq, wt_seq) if a == b)
        if match > best_match:
            best_match = match
            best_wt = wt_name
    n_mut = len(seq) - best_match if best_wt else "?"
    pct_id = best_match / len(seq) * 100 if best_wt else 0
    print(f"Seq {row['Seq_ID']}: closest={best_wt}, muts={n_mut}, identity={pct_id:.1f}%")
