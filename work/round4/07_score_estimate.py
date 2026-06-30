"""
Round 4 比赛得分估算 - 基于文献先验 + ESMFold 数据
==================================
得分规则:
  - 相对亮度 = Finitial / Finitial_WT  (WT = sfGFP, baseline=1.0)
  - 热稳定性 = Ffinal / Finitial  (72°C 后保留比例)
  - 综合得分 = 相对亮度 × 热稳定性
  - Finitial < 0.3 × Finitial_WT 直接淘汰为 0 分

文献先验:
  - sfGFP WT Tm: 86.1°C (CD denaturation, BMC Res Notes 2023)
  - sfGFP 72°C 半小时保留率: ~85% (Tm-T=14°C 余量大)
  - mBaoJin Tm: ~92°C, 70°C 加热 1h 完全保持
  - StayGold/mBaoJin 在 72°C 应几乎 100% 保留
  - I152S/Q69L 等点突变对 brightness 小幅扰动 (±20%)
  - htFuncLib sf:acid.3 报告 Tm 可达 96°C

热失活动力学:
  - 1阶 Arrhenius 模型: k(T) = k_ref * exp(Ea/R * (1/T_ref - 1/T))
  - sfGFP 在 90°C 半衰期 ~30 min
  - 假设 72°C 加热时长 30 min (官方未给, 假设值)
  - Ffinal/Finit = exp(-k*t)
"""
import json, math
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
OUT = ROOT / "work" / "round4"

with open(OUT / "final_6_round4_v2.json", encoding="utf-8") as f:
    final_6 = json.load(f)

# ============================================================
# 文献先验参数
# ============================================================
# sfGFP 作为 WT 参考
SFGFP_WT_BRIGHTNESS = 1.0  # 定义为 1.0
SFGFP_WT_TM = 86.1  # °C (BMC Res Notes 2023)
TEST_T = 72  # 比赛热处理温度
HEAT_TIME_MIN = 30  # 假设热处理 30 分钟 (官方未明确)

# 突变效应先验 (来自论文/数据集)
mutation_effects = {
    # sfGFP 已含 (基线参考)
    # 单点突变
    "I152S": {"brightness_factor": 1.05, "dTm": +2},  # Round 1 验证, chromophore 邻位优化
    "Q69L": {"brightness_factor": 1.10, "dTm": +6},   # htFuncLib sf:acid 验证
    "S72A": {"brightness_factor": 1.05, "dTm": +4},   # htFuncLib
    "T108V": {"brightness_factor": 1.00, "dTm": +3},  # htFuncLib
    "K166V": {"brightness_factor": 0.95, "dTm": +1},  # cgre-style hotspot
    "S30R": {"brightness_factor": 1.00, "dTm": +5},   # sfGFP 已含, 在avGFP+8°C
    "F64L": {"brightness_factor": 1.20, "dTm": +3},   # sfGFP 核心
    "F99S": {"brightness_factor": 1.15, "dTm": +4},
    "M153T": {"brightness_factor": 1.10, "dTm": +3},
    "V163A": {"brightness_factor": 1.05, "dTm": +2},
    "S65T": {"brightness_factor": 1.30, "dTm": +2},   # chromophore 成熟加速
}


def estimate_finit_relative(c):
    """
    估算 Finitial / Finitial_sfGFP_WT
    基于:
      1. 骨架 baseline brightness
      2. 突变累加效应
      3. pLDDT 折叠概率 (惩罚)
    """
    # Step 1: 骨架基础亮度 (相对 sfGFP WT)
    scaffold_baseline = {
        "sfGFP": 1.00,        # 定义为 1.0 WT
        "avGFP": 0.90,        # 略低于 sfGFP (折叠差), 但 sfGFP 突变叠加可超过
        "amacGFP": 0.85,      # baseline 3.97 vs sfGFP 4.0, 加 sfGFP 突变可平
        "cgreGFP": 1.10,      # baseline 4.50 高
        "mBaoJin": 1.20,      # 报告比 mNeonGreen 亮 (~mNG 1.3-1.4× sfGFP)
    }
    base = scaffold_baseline.get(c["scaffold"], 1.0)

    # Step 2: 突变累加 (从 notes 提取已知突变效应)
    notes = c.get("notes", "")
    mut_factor = 1.0
    # sfGFP 风格突变 (S65T/F64L/F99S/M153T/V163A 等)
    if "sfGFP" in c["name"] and c["scaffold"] != "sfGFP":
        # avGFP/amacGFP + sfGFP 风格 => 折叠改善
        mut_factor *= 1.5  # 大幅提升
    
    # I152S 加成 (Round 1 验证)
    if "I152S" in notes or "I152S" in c["name"]:
        mut_factor *= 1.05
    
    # htFuncLib Q69L 加成
    if "Q69L" in c["name"] or "acid" in c["name"].lower():
        mut_factor *= 1.10
    
    # mBaoJin 单体表面突变 (微小影响)
    if c["scaffold"] == "mBaoJin":
        mut_factor *= 0.98  # 略有风险

    # Step 3: pLDDT 折叠概率惩罚
    # ESMFold pLDDT 40-50 视为"中等置信", 几乎都能折叠 (实证GFP系列)
    # pLDDT 35 以下风险高, 30 以下大概率失败
    plddt = c["plddt_mean"]
    if plddt >= 45:
        fold_prob = 0.95
    elif plddt >= 40:
        fold_prob = 0.85
    elif plddt >= 35:
        fold_prob = 0.65  # 不确定 (mBaoJin 在此区间)
    else:
        fold_prob = 0.30  # 高风险

    # chromophore 区域 pLDDT 单独惩罚 (更直接反映 chromophore 形成)
    chromo_plddt = c["plddt_chromo_region"]
    if chromo_plddt < 35:
        chromo_factor = 0.70  # mBaoJin cb=35 风险
    elif chromo_plddt < 45:
        chromo_factor = 0.90
    else:
        chromo_factor = 1.00

    finit_rel = base * mut_factor * fold_prob * chromo_factor
    return finit_rel, {
        "scaffold_base": base,
        "mut_factor": round(mut_factor, 2),
        "fold_prob": fold_prob,
        "chromo_factor": chromo_factor,
    }


def estimate_thermal_retention(c, heat_time_min=30, test_T=72):
    """
    估算 72°C 加热后的保留率 Ffinal / Finit
    使用 1阶动力学: ln(F/F0) = -k * t
    k(T) 由 Tm 决定 (T < Tm 时极慢, T > Tm 时快速)
    """
    tm = c["expected_tm"]

    # 根据 Tm 与 test_T 的关系估算保留率
    # T < Tm-15: 几乎完全保留 (>95%)
    # Tm-15 < T < Tm-5: 中度保留 (70-95%)
    # Tm-5 < T < Tm: 显著损失 (30-70%)
    # T > Tm: 大量损失 (<30%)
    delta = tm - test_T  # 余量

    if delta >= 20:
        retention = 0.98  # 几乎完全 (mBaoJin Tm=92, 余量20°C)
    elif delta >= 14:
        retention = 0.92  # sfGFP Tm=86, 余量14°C
    elif delta >= 10:
        retention = 0.85
    elif delta >= 6:
        retention = 0.70
    elif delta >= 0:
        retention = 0.50
    else:
        retention = 0.20

    return retention, {"delta_Tm_T": delta}


# ============================================================
# 估算所有候选
# ============================================================
print("=" * 110)
print(f"{'Seq':<4} {'name':<32} {'scaffold':<8} {'Tm':>3} {'pLDDT':>5} "
      f"{'Finit_rel':>9} {'Therm_ret':>9} {'综合分':>6} {'是否>0.3':>8}")
print("=" * 110)

estimates = []
for c in final_6:
    finit_rel, finit_detail = estimate_finit_relative(c)
    therm_ret, therm_detail = estimate_thermal_retention(c)
    
    # 综合分
    combined = finit_rel * therm_ret
    pass_threshold = finit_rel >= 0.3
    
    if not pass_threshold:
        combined = 0  # 淘汰

    icon = "✅" if pass_threshold else "❌"
    print(f"{c['Seq_ID']:<4} {c['name']:<32} {c['scaffold']:<8} "
          f"{c['expected_tm']:>3} {c['plddt_mean']:>5.1f} "
          f"{finit_rel:>9.2f} {therm_ret:>9.2f} {combined:>6.2f} {icon:>8}")

    estimates.append({
        "Seq_ID": c["Seq_ID"],
        "name": c["name"],
        "scaffold": c["scaffold"],
        "expected_tm": c["expected_tm"],
        "plddt_mean": c["plddt_mean"],
        "finit_relative": round(finit_rel, 3),
        "thermal_retention": round(therm_ret, 3),
        "combined_score": round(combined, 3),
        "pass_threshold": pass_threshold,
        "details": {**finit_detail, **therm_detail},
    })

# 找最高分 (Best Top-1)
top1 = max(estimates, key=lambda x: x["combined_score"])
print("\n" + "=" * 80)
print(f"🏆 预测 Best Top-1: Seq {top1['Seq_ID']} - {top1['name']}")
print(f"   综合分: {top1['combined_score']}")
print(f"   = Finit_rel {top1['finit_relative']} × Thermal_ret {top1['thermal_retention']}")
print("=" * 80)

# 多场景分析 (乐观/悲观)
print("\n" + "=" * 80)
print("场景分析 (考虑不确定性):")
print("=" * 80)

# 悲观: 所有 mut_factor 打 0.7 折, pLDDT 阈值更严
def pessimistic(c):
    finit_rel, _ = estimate_finit_relative(c)
    therm_ret, _ = estimate_thermal_retention(c)
    return finit_rel * 0.65 * therm_ret * 0.85  # 悲观折扣

# 乐观: 所有 mut_factor 涨 1.3, 热稳更好
def optimistic(c):
    finit_rel, _ = estimate_finit_relative(c)
    therm_ret, _ = estimate_thermal_retention(c)
    return finit_rel * 1.4 * therm_ret * 1.10  # 乐观加成

print(f"{'Seq':<4} {'name':<32} {'悲观':>6} {'中性':>6} {'乐观':>6}")
print("-" * 65)
for c, est in zip(final_6, estimates):
    pess = pessimistic(c)
    neu = est["combined_score"]
    opt = optimistic(c)
    if not est["pass_threshold"]:
        pess = neu = opt = 0
    print(f"{c['Seq_ID']:<4} {c['name']:<32} {pess:>6.2f} {neu:>6.2f} {opt:>6.2f}")

# Best Top-1 三个场景
print("\n" + "=" * 80)
print("🎯 Best Top-1 三场景预测:")
print("=" * 80)
pess_scores = [pessimistic(c) if est["pass_threshold"] else 0 for c, est in zip(final_6, estimates)]
neu_scores = [est["combined_score"] for est in estimates]
opt_scores = [optimistic(c) if est["pass_threshold"] else 0 for c, est in zip(final_6, estimates)]

print(f"  🔴 悲观 (折叠较差/亮度低估): Best Top-1 = {max(pess_scores):.2f}")
print(f"  🟡 中性 (基于文献先验):     Best Top-1 = {max(neu_scores):.2f}")
print(f"  🟢 乐观 (mut 协同 + 高保留): Best Top-1 = {max(opt_scores):.2f}")

# 保存
with open(OUT / "score_estimates.json", "w", encoding="utf-8") as f:
    json.dump({
        "estimates": estimates,
        "best_top1": top1,
        "scenarios": {
            "pessimistic_top1": round(max(pess_scores), 3),
            "neutral_top1": round(max(neu_scores), 3),
            "optimistic_top1": round(max(opt_scores), 3),
        },
        "assumptions": {
            "sfGFP_WT_baseline": 1.0,
            "sfGFP_WT_Tm": 86.1,
            "test_temperature": 72,
            "heat_time_min": 30,
            "fold_prob_from_plddt": "45+ → 0.95; 40-45 → 0.85; 35-40 → 0.65; <35 → 0.30",
        }
    }, f, indent=2, ensure_ascii=False)

print(f"\n✓ 保存到 work/round4/score_estimates.json")

print("\n" + "=" * 80)
print("⚠️ 重要免责声明:")
print("=" * 80)
print("""
1. 这是基于文献先验的估算, 不是实验值
2. 关键不确定性:
   - 实际比赛 CFPS 体系条件未知
   - 72°C 加热时长未知 (假设 30min)
   - 突变协同效应 (epistasis) 难以预测
   - mBaoJin 在 CFPS 中表达效率未知
3. 估算误差范围可能 ±50%
4. 实际比赛得分需等官方测定
""")
