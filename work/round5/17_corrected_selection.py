"""
Round 5 修正版: 以结构完整性 (pTM/pLDDT) 为主排序
==================================
用户指出的核心问题:
  1. mBaoJin pTM~0.36 = 无序卷曲, 不是蛋白质
  2. Finit/Therm 是手工估值, 与 pLDDT/pTM 矛盾
  3. 手工点突变破坏上位协同网络
  4. 比 Round 4 (pTM 0.765) 倒退
  5. 排序优先级搞反了

正确策略:
  - pTM > 0.5 是"可能正确折叠"的最低门槛
  - pLDDT > 45 是第二门槛
  - 在结构完整的前提下, 突变数越少越好
  - Tm 文献先验仅作辅助参考
"""
import json, pandas as pd
from pathlib import Path
from collections import defaultdict

ROOT = Path(r"D:\生信\2026Protein Design")
R5 = ROOT / "work" / "round5"

# 加载所有候选 (104 条, 含 EVOLVEpro)
with open(R5 / "evolvepro_scored.json", encoding="utf-8") as f:
    all_c = json.load(f)

print(f"总候选: {len(all_c)}")

# ============================================================
# 正确的筛选逻辑
# ============================================================
# 第一关: 结构完整性 (pTM >= 0.5 AND pLDDT >= 45)
# 这是"可能正确折叠"的最低门槛
print("\n" + "=" * 100)
print("第一关: 结构完整性筛选 (pTM >= 0.5 AND pLDDT >= 45)")
print("=" * 100)

struct_ok = [r for r in all_c if (r["ptm"] or 0) >= 0.5 and r["plddt_mean"] >= 45]
print(f"通过: {len(struct_ok)} / {len(all_c)}")

# 第二关: 30% 阈值安全 (突变数 ≤ 12)
# 突变数 ≤ 12 的候选大概率保留 >30% WT 亮度
print(f"\n第二关: 30% 阈值安全 (突变数 ≤ 12)")
safe = [r for r in struct_ok if r["n_muts"] <= 12]
print(f"通过: {len(safe)} / {len(struct_ok)}")

# 如果安全候选不够, 放宽到 mut ≤ 15
if len(safe) < 6:
    safe_15 = [r for r in struct_ok if r["n_muts"] <= 15]
    print(f"放宽到 mut ≤ 15: {len(safe_15)}")

# 第三关: 排序 — 以 pTM 为主 (结构完整性), 突变数为辅
safe.sort(key=lambda x: (-(x["ptm"] or 0), x["n_muts"]))

print(f"\n{'='*120}")
print(f"结构完整 + 安全候选 (pTM≥0.5, pLDDT≥45, mut≤12), 按 pTM 降序:")
print(f"{'='*120}")
print(f"{'#':<3} {'name':<35} {'scaffold':<14} {'mut':>3} {'pLDDT':>5} {'cb':>5} {'pTM':>5} {'evo':>5}")
print("-" * 110)
for i, r in enumerate(safe[:25], 1):
    print(f"{i:<3} {r['name'][:35]:<35} {r['scaffold'][:14]:<14} {r['n_muts']:>3} "
          f"{r['plddt_mean']:>5.1f} {r['plddt_chromo_region']:>5.1f} {r['ptm'] or 0:>5.3f} "
          f"{r.get('evolvepro_pred',0):>5.2f}")

# ============================================================
# 也看 MPNN/LMPNN 高 pTM 候选 (即使突变多)
# ============================================================
print(f"\n{'='*120}")
print(f"高 pTM MPNN/LMPNN 候选 (pTM≥0.6, 突变数不限):")
print(f"{'='*120}")
mpnn_high = [r for r in all_c if (r["ptm"] or 0) >= 0.6 and ("MPNN" in r["scaffold"] or "LMPNN" in r["scaffold"])]
mpnn_high.sort(key=lambda x: -(x["ptm"] or 0))
print(f"{'#':<3} {'name':<35} {'scaffold':<14} {'mut':>3} {'pLDDT':>5} {'pTM':>5}")
print("-" * 80)
for i, r in enumerate(mpnn_high[:15], 1):
    print(f"{i:<3} {r['name'][:35]:<35} {r['scaffold'][:14]:<14} {r['n_muts']:>3} "
          f"{r['plddt_mean']:>5.1f} {r['ptm'] or 0:>5.3f}")

# ============================================================
# 混合策略 Top-6: 2 条 MPNN 冲刺 + 4 条安全保底
# ============================================================
print(f"\n{'='*120}")
print(f"🎯 混合策略 Top-6 (2 MPNN 冲刺 + 4 安全保底)")
print(f"{'='*120}")

# 冲刺: MPNN/LMPNN 高 pTM (pTM >= 0.65)
sprint = [r for r in mpnn_high if (r["ptm"] or 0) >= 0.65][:2]
# 保底: 安全候选 (pTM >= 0.5, mut <= 12), 不同骨架
# 排除已选冲刺的骨架
sprint_scaffolds = {r["scaffold"] for r in sprint}
保底_pool = [r for r in safe if r["scaffold"] not in sprint_scaffolds]

# 多样性选择保底
priority_safe = ["sfGFP", "avGFP", "amacGFP", "mBaoJin", "cgreGFP"]
保底_selected = []
seen = set()
sc_cnt = defaultdict(int)
for sc in priority_safe:
    for r in 保底_pool:
        if r["scaffold"] == sc and r["seq"] not in seen:
            保底_selected.append(r); seen.add(r["seq"]); sc_cnt[sc] += 1; break
# 补足到 4
for r in 保底_pool:
    if len(保底_selected) >= 4: break
    if r["seq"] in seen: continue
    if sc_cnt[r["scaffold"]] >= 2: continue
    保底_selected.append(r); seen.add(r["seq"]); sc_cnt[r["scaffold"]] += 1

final_6 = sprint[:2] + 保底_selected[:4]

# 按 pTM 排序 (结构完整性为主)
final_6.sort(key=lambda x: -(x["ptm"] or 0))

print(f"\n{'Seq':<4} {'name':<35} {'scaffold':<14} {'mut':>3} {'pLDDT':>5} {'cb':>5} {'pTM':>5} {'角色'}")
print("-" * 100)
for i, r in enumerate(final_6, 1):
    role = "冲刺" if r in sprint else "保底"
    print(f"{i:<4} {r['name'][:35]:<35} {r['scaffold'][:14]:<14} {r['n_muts']:>3} "
          f"{r['plddt_mean']:>5.1f} {r['plddt_chromo_region']:>5.1f} {r['ptm'] or 0:>5.3f} {role}")

# 保存
output = []
for i, r in enumerate(final_6, 1):
    output.append({
        "Seq_ID": i, "name": r["name"], "scaffold": r["scaffold"],
        "n_muts": r["n_muts"], "expected_tm": r.get("expected_tm", 80),
        "plddt_mean": r["plddt_mean"], "plddt_chromo_region": r["plddt_chromo_region"],
        "ptm": r["ptm"], "evolvepro_pred": r.get("evolvepro_pred"),
        "notes": r.get("notes",""), "seq": r["seq"],
    })
with open(R5 / "final_6_round5_hybrid.json", "w", encoding="utf-8") as f:
    json.dump(output, f, indent=2, ensure_ascii=False)

sub = pd.DataFrame({
    "Team_Name": ["YourTeamName"] * 6,
    "Seq_ID": list(range(1, 7)),
    "Sequence": [r["seq"] for r in final_6]
})
sub_path = R5 / "submission_round5_hybrid.csv"
sub.to_csv(sub_path, index=False)

# 合规
excl = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())
all_ok = all(s.startswith("M") and 220<=len(s)<=250 and not (set(s)-set("ACDEFGHIKLMNPQRSTVWY")) and s not in excl for s in sub["Sequence"])
print(f"\n📋 合规: {'✅' if all_ok else '❌'}")
print(f"📁 提交: {sub_path}")

# 各候选分析
print(f"\n{'='*100}")
print("📊 各候选分析:")
print(f"{'='*100}")
for i, r in enumerate(final_6, 1):
    plddt = r["plddt_mean"]; ptm = r["ptm"] or 0; n_mut = r["n_muts"]
    role = "冲刺" if r in sprint else "保底"
    
    # 结构完整性
    if ptm >= 0.7:
        struct = "✅ 高置信折叠"
    elif ptm >= 0.5:
        struct = "🟡 可能正确折叠"
    else:
        struct = "🔴 折叠不确定"
    
    # 30% 阈值
    if n_mut <= 5:
        threshold = "✅ 安全 (>80% 亮度保留)"
    elif n_mut <= 12:
        threshold = "🟡 中等 (~60% 亮度保留)"
    elif n_mut <= 60:
        threshold = "🟠 风险 (~30% 亮度保留, 阈值边缘)"
    else:
        threshold = "🔴 高风险 (<30% 亮度, 可能淘汰)"
    
    print(f"\n  Seq {i}: {r['name'][:35]} ({role})")
    print(f"    结构: pLDDT={plddt:.1f} pTM={ptm:.3f} → {struct}")
    print(f"    亮度: mut={n_mut} → {threshold}")
    print(f"    热稳: Tm~{r.get('expected_tm',80)}°C")
