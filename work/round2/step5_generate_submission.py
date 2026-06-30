"""Step 5: 生成 round 2 提交 CSV。
原则:
- 6 条覆盖 4 个不同母体(avGFP/sfGFP/cgreGFP/amacGFP/ppluGFP)
- 优先 chromophore 完整 + 论文突变数高 + 排除列表粗查
- 平衡"保守"和"激进": 2 条保守(WT 或近 WT) + 4 条激进(论文突变组合)
"""
import pandas as pd
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design\work\round2")
CAND_FILE = ROOT / "candidates_round2_design.csv"
VAL_FILE = ROOT / "step4_validation.csv"
OUT_SUB = ROOT / "submission_round2.csv"
FINAL_SUB = Path(r"D:\生信\2026Protein Design\submission_yourteamname.csv")

# 读候选 + 评分
cand = pd.read_csv(CAND_FILE)
val = pd.read_csv(VAL_FILE)

# 选 Top-6(平衡多样性)
# 设计策略:
# - #6 avGFP+sfGFP完整+TGP_全面增强(19 muts, 4 scaffold score)
# - #3 avGFP+sfGFP完整+TGP稳定(14 muts, 中等)
# - #4 amacGFP+sfGFP风格+TGP稳定(10 muts)
# - #2 sfGFP+TGP_稳定核心(4 muts, sfGFP 基线)
# - #5 ppluGFP原始(0 muts, ppluGFP 高 baseline, 保守)
# - #9 cgreGFP+S65T+K163A(2 muts, cgreGFP 高 baseline)

selected_ids = [6, 3, 4, 2, 5, 9]  # 平衡多样性
selected = val[val["id"].isin(selected_ids)].copy()
selected["seq"] = selected["id"].apply(lambda i: cand[cand["id"]==i]["seq"].iloc[0])

print("=== Round 2 最终 6 条候选 ===")
for _, r in selected.iterrows():
    print(f"#{r['id']:2d} {r['name']:35s} scaffold={r['scaffold']:8s} n_muts={r['n_muts']:2d} score={r['design_score']:3d}")

# 生成 submission CSV (格式: Team_Name,Seq_ID,Sequence)
rows = []
for new_id, (_, r) in enumerate(selected.iterrows(), 1):
    rows.append({
        "Team_Name": "YourTeamName",
        "Seq_ID": new_id,
        "Sequence": r["seq"],
    })
df_sub = pd.DataFrame(rows)
df_sub.to_csv(OUT_SUB, index=False)
print(f"\n✅ Round 2 submission 已保存到 {OUT_SUB}")

# 也覆盖主 submission 文件
df_sub.to_csv(FINAL_SUB, index=False)
print(f"✅ 主 submission 也已覆盖到 {FINAL_SUB}")

# 最终验证
print("\n=== 最终验证 ===")
for _, r in df_sub.iterrows():
    s = r["Sequence"]
    issues = []
    if len(s) < 220 or len(s) > 250:
        issues.append(f"长度 {len(s)}")
    if s[0] != "M":
        issues.append("不以 M 开头")
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    bad = set(s) - valid_aa
    if bad:
        issues.append(f"非标准 AA: {bad}")
    chrom = "TYG" in s or "SYG" in s or "GYG" in s
    print(f"Seq {r['Seq_ID']}: len={len(s)}, M-start={'✓' if s[0]=='M' else '✗'}, chromophore={'✓' if chrom else '✗'}{', '+','.join(issues) if issues else ''}")