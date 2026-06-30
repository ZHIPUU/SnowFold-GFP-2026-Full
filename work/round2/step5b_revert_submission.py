"""回退到 paper-knowledge-based Top-6(epistasis 模型对 OOD 候选不可靠)"""
import pandas as pd
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design\work\round2")
CAND_FILE = ROOT / "candidates_round2_design.csv"
OUT_SUB = Path(r"D:\生信\2026Protein Design\submission_yourteamname.csv")

cand = pd.read_csv(CAND_FILE)

# Paper-knowledge Top-6(基于论文先验 + scaffold 多样性 + 设计合理性)
selected_ids = [6, 3, 4, 2, 5, 9]  # 与 step5 一致

rows = []
for new_id, cid in enumerate(selected_ids, 1):
    row = cand[cand["id"] == cid].iloc[0]
    rows.append({
        "Team_Name": "YourTeamName",
        "Seq_ID": new_id,
        "Sequence": row["seq"],
    })

df_sub = pd.DataFrame(rows)
df_sub.to_csv(OUT_SUB, index=False)
print(f"✅ Submission 已回退到 paper-knowledge Top-6 → {OUT_SUB}")

# 验证
print("\n=== 最终提交 ===")
for new_id, cid in enumerate(selected_ids, 1):
    row = cand[cand["id"] == cid].iloc[0]
    s = row["seq"]
    chrom = "TYG" in s or "SYG" in s or "GYG" in s
    print(f"  Seq {new_id}: {row['name']:35s} scaffold={row['scaffold']:8s} len={len(s)} chromo={'✓' if chrom else '✗'}")