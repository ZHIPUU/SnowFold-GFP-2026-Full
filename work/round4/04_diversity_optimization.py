"""
Round 4 优化 Step 1: 解决多样性问题 + 扩展骨架
==================================
问题诊断 (来自 diagnose_diversity.py):
  - 6条提交中 5条几乎相同 (汉明距离 1-4)
  - sfGFP 占 5/6, 系统性风险
  - mBaoJin 仅 1 条, 失去热稳奖竞争力

优化策略:
  - 增加 mBaoJin/cgreGFP/avGFP 多骨架候选
  - 多机制突变组合 (htFuncLib + S30R + I152S + 表面优化)
  - 确保最终 6 条来自至少 3 种骨架
"""
import re, json, pandas as pd
from pathlib import Path
from collections import Counter

ROOT = Path(r"D:\生信\2026Protein Design")
OUT = ROOT / "work" / "round4"

with open(ROOT / "AAseqs of 5 GFP proteins_20260511.txt") as f:
    wt_text = f.read()
wts = {}
for block in wt_text.split(">"):
    block = block.strip()
    if not block: continue
    lines = block.split("\n")
    name = lines[0].strip()
    seq = "".join(l.strip() for l in lines[1:] if l.strip() and not l.startswith("#"))
    wts[name] = seq

avGFP = wts["avGFP"]
sfGFP = wts["sfGFP"]
amacGFP = wts["amacGFP"]
cgreGFP = wts["cgreGFP"]
mbaojin = "MVSKGEEENMASTPFKFQLKGTINGKSFTVEGEGEGNSHEGSHKGKYVCTSGKLPMSWAALGTTFGYGMKYYTKYPSGLKNWFREVMPGGFTYDRHIQYKGDGSIHAKHQHFMKNGTYHNIVEFTGQDFKENSPVLTGDMNVSLPNEVPQIPRDDGVECPVTLLYPLLSDKSKYVEAHQYTICKPLHNQPAPDVPYHWIRKQYTQSKDDAEERDHICQSETLEAHLKGMDELYK"

def apply(seq, muts):
    s = list(seq)
    for m in muts:
        match = re.match(r"([A-Z])(\d+)([A-Z])", m)
        if not match: continue
        fr, pos, to = match.group(1), int(match.group(2)), match.group(3)
        idx = pos - 1
        if idx >= len(s): continue
        if s[idx] != fr:
            # 严格模式: 不应用错误突变
            continue
        s[idx] = to
    return "".join(s)

def n_muts(seq, parent):
    return sum(1 for a, b in zip(seq, parent) if a != b)

def chromo(seq):
    for t in ["TYG", "SYG", "GYG", "CYG", "HYG"]:
        if t in seq: return t
    return "?"

# ============================================================
# 新候选 (多骨架 + 多机制)
# ============================================================
candidates = []

# --------------------------------------------------------------------
# X1-X5: avGFP 系列 (完整不同变体, 避免与 sfGFP 重复)
# --------------------------------------------------------------------
# avGFP pos 检查: 30=S, 64=L(已sfGFP化), 65=S, 99=F, 152=I, 153=M, 163=V, 65=S
# avGFP + sfGFP核心4 (不含 F64L 因已有) = S65T, F99S, M153T, V163A
# avGFP + 完整 sfGFP (除 F64L)
# avGFP + sfGFP10 + I152S + Q69L (Round 1 9.50 + htFuncLib)
sfgfp_5core_av = ["S65T", "F99S", "M153T", "V163A"]
sfgfp_full_av = sfgfp_5core_av + ["S30R", "Y39N", "N105T", "Y145F", "I171V", "A206V"]

# X1: avGFP + sfGFP 4 核心 (最保守起点)
candidates.append({
    "name": "X1_avGFP_sf4core",
    "seq": apply(avGFP, sfgfp_5core_av),
    "parent": avGFP, "scaffold": "avGFP",
    "role": "diversity_avGFP", "expected_tm": 78,
    "notes": "avGFP + sfGFP 4 折叠核心, 4mut, 最保守"
})

# X2: avGFP + sfGFP 完整 + I152S (Round 1 Seq 6 重现, 综合分 9.50 估计)
candidates.append({
    "name": "X2_avGFP_sfGFP_I152S",
    "seq": apply(avGFP, sfgfp_full_av + ["I152S"]),
    "parent": avGFP, "scaffold": "avGFP",
    "role": "diversity_avGFP", "expected_tm": 82,
    "notes": "avGFP + sfGFP10 + I152S (Round 1 Seq 6 估 9.50), 11mut"
})

# X3: avGFP + sfGFP10 + Q69L + S72A (htFuncLib 双突变叠加)
candidates.append({
    "name": "X3_avGFP_sfGFP_acid2",
    "seq": apply(avGFP, sfgfp_full_av + ["Q69L", "S72A"]),
    "parent": avGFP, "scaffold": "avGFP",
    "role": "diversity_avGFP", "expected_tm": 88,
    "notes": "avGFP + sfGFP10 + Q69L + S72A htFuncLib, 12mut"
})

# X4: avGFP + sfGFP10 + I152S + Q69L (综合)
candidates.append({
    "name": "X4_avGFP_sfGFP_I152S_Q69L",
    "seq": apply(avGFP, sfgfp_full_av + ["I152S", "Q69L"]),
    "parent": avGFP, "scaffold": "avGFP",
    "role": "diversity_avGFP", "expected_tm": 90,
    "notes": "avGFP + sfGFP10 + I152S + Q69L 三机制叠加, 12mut"
})

# --------------------------------------------------------------------
# Y1-Y3: cgreGFP 系列 (新骨架, Round 3 试过但 pLDDT 低)
# cgreGFP baseline brightness 最高 (4.50), 但需谨慎突变
# Round 1 验证 cgreGFP + K163S + G130C + A8Q + W192A (Seq 3)
# --------------------------------------------------------------------
# cgreGFP pos 检查
print(f"cgreGFP pos: 8={cgreGFP[7]}, 130={cgreGFP[129]}, 163={cgreGFP[162]}, 192={cgreGFP[191]}, 65={cgreGFP[64]}")

# Y1: cgreGFP + S65T (chromophore 现代化)
candidates.append({
    "name": "Y1_cgreGFP_S65T",
    "seq": apply(cgreGFP, ["S65T"]),
    "parent": cgreGFP, "scaffold": "cgreGFP",
    "role": "diversity_cgreGFP", "expected_tm": 78,
    "notes": "cgreGFP + S65T 单点 chromophore 现代化, baseline 4.50 高, 1mut"
})

# Y2: cgreGFP + S65T + K163A (Round 3 c gre 候选)
candidates.append({
    "name": "Y2_cgreGFP_S65T_K163A",
    "seq": apply(cgreGFP, ["S65T", "K163A"]),
    "parent": cgreGFP, "scaffold": "cgreGFP",
    "role": "diversity_cgreGFP", "expected_tm": 80,
    "notes": "cgreGFP + S65T + K163A, 2mut"
})

# Y3: cgreGFP + Round1验证组合 (K163S/G130C/A8Q/W192A)
candidates.append({
    "name": "Y3_cgreGFP_R1_combo",
    "seq": apply(cgreGFP, ["K163S", "G130C", "A8Q", "W192A"]),
    "parent": cgreGFP, "scaffold": "cgreGFP",
    "role": "diversity_cgreGFP", "expected_tm": 80,
    "notes": "cgreGFP + Round1 Seq 3 组合, 4mut, baseline 4.50"
})

# --------------------------------------------------------------------
# Z1-Z3: amacGFP 系列 (跨骨架多样性)
# amacGFP baseline 3.97, K166I hotspot
# --------------------------------------------------------------------
print(f"amacGFP pos: 65={amacGFP[64]}, 166={amacGFP[165]}, 65 chromophore={amacGFP[64:67]}")

# Z1: amacGFP + sfGFP 风格 (Round 3 验证候选)
candidates.append({
    "name": "Z1_amacGFP_sfGFP5",
    "seq": apply(amacGFP, ["S65T", "F99S", "M153T", "I166T", "I171V"]),
    "parent": amacGFP, "scaffold": "amacGFP",
    "role": "diversity_amacGFP", "expected_tm": 78,
    "notes": "amacGFP + sfGFP 风格 5mut (含 amacGFP hotspot K166I)"
})

# Z2: amacGFP + K166I (单点 hotspot)
candidates.append({
    "name": "Z2_amacGFP_K166I",
    "seq": apply(amacGFP, ["K166I"]),
    "parent": amacGFP, "scaffold": "amacGFP",
    "role": "diversity_amacGFP", "expected_tm": 75,
    "notes": "amacGFP + K166I hotspot 单突变, 1mut"
})

# Z3: amacGFP + sfGFP 风格 + I152S
candidates.append({
    "name": "Z3_amacGFP_sfGFP5_I152S",
    "seq": apply(amacGFP, ["S65T", "F99S", "M153T", "I166T", "I171V", "I152S"]),
    "parent": amacGFP, "scaffold": "amacGFP",
    "role": "diversity_amacGFP", "expected_tm": 80,
    "notes": "amacGFP + sfGFP5 + I152S, 6mut"
})

# --------------------------------------------------------------------
# P1-P5: mBaoJin 系列扩展 (5 条不同, 选最佳1-2)
# --------------------------------------------------------------------
# 已有: M1-M8 已验证, 选最佳 K173R / Y196F / 组合
# 新增: 多位点组合
candidates.append({
    "name": "P1_mBaoJin_K173R_V194I_Y196F",
    "seq": apply(mbaojin, ["K173R", "V194I", "Y196F"]),
    "parent": mbaojin, "scaffold": "mBaoJin",
    "role": "thermostable_hero", "expected_tm": 92,
    "notes": "mBaoJin + 3 表面保守突变 (K173R/V194I/Y196F), 3mut"
})

candidates.append({
    "name": "P2_mBaoJin_K173R_L222M",
    "seq": apply(mbaojin, ["K173R", "L222M"]),
    "parent": mbaojin, "scaffold": "mBaoJin",
    "role": "thermostable_hero", "expected_tm": 92,
    "notes": "mBaoJin + K173R + L222M C端保守, 2mut"
})

candidates.append({
    "name": "P3_mBaoJin_V194I_D230E",
    "seq": apply(mbaojin, ["V194I", "D230E"]),
    "parent": mbaojin, "scaffold": "mBaoJin",
    "role": "thermostable_hero", "expected_tm": 92,
    "notes": "mBaoJin + V194I + D230E 双保守, 2mut"
})

# --------------------------------------------------------------------
# 验证 + 去重 + 加载已有 22 条
# --------------------------------------------------------------------
excl_seqs = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())

print("\n" + "=" * 90)
print(f"{'name':<32} {'scaffold':<10} {'mut':>3} {'chromo':>6} {'len':>4} {'excl':>5} {'status'}")
print("=" * 90)

valid_new = []
for c in candidates:
    seq = c["seq"]
    nm = n_muts(seq, c["parent"])
    ch = chromo(seq)
    in_excl = seq in excl_seqs
    issues = []
    if not seq.startswith("M"): issues.append("no-M")
    if not (220 <= len(seq) <= 250): issues.append(f"len={len(seq)}")
    if ch == "?": issues.append("no-chromo")
    if nm == 0: issues.append("zero-mut")

    excl_mark = "❌" if in_excl else "✓"
    status = "✓" if not issues and not in_excl else ",".join(issues) + (" EXCL!" if in_excl else "")
    print(f"{c['name']:<32} {c['scaffold']:<10} {nm:>3} {ch:>6} {len(seq):>4} {excl_mark:>5} {status}")

    c["length"] = len(seq)
    c["n_muts"] = nm
    c["chromophore"] = ch
    c["in_exclusion"] = in_excl
    c["valid"] = (not issues) and (not in_excl)
    if c["valid"]:
        valid_new.append(c)

print(f"\n新增有效: {len(valid_new)} / {len(candidates)}")

# 合并扩展候选池
with open(OUT / "candidates_round4_extended.json", encoding="utf-8") as f:
    existing = json.load(f)

all_seqs = set(c["seq"] for c in existing)
added = 0
for c in valid_new:
    if c["seq"] in all_seqs: continue
    c_save = {k: v for k, v in c.items() if k != "parent"}
    c_save["parent_scaffold"] = c["scaffold"]
    existing.append(c_save)
    all_seqs.add(c["seq"])
    added += 1

print(f"\n合并: 旧 {len(existing) - added} + 新增 {added} = 总 {len(existing)}")

with open(OUT / "candidates_round4_v2.json", "w", encoding="utf-8") as f:
    json.dump(existing, f, indent=2, ensure_ascii=False)

print("\n保存到 candidates_round4_v2.json")
print("\n按骨架统计:")
for sc, n in Counter(c.get("scaffold") for c in existing).most_common():
    print(f"  {sc}: {n} 条")
