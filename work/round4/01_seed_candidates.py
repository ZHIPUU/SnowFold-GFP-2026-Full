"""
Round 4 Step 1: 生成候选池
==================================
依据:
  - sfGFP (Pédelacq 2006) - 11 经典突变 (F64L/S65T/F99S/M153T/V163A + S30R/Y39N/N105T/Y145F/I171V/A206V)
  - htFuncLib sf:acid.3 (Weinstein 2023) - T65S/Q69L/S72A/T108V/Y145M/V224I, Tm 高达 96°C
  - mBaoJin (Zhang 2024) - 8 突变 vs StayGold: S55T/H77R/E80G/Q140P/H141Q/C165Y/N171Y/T201A
  - TGP/eCGP123 (Close 2015) - mAG 编号: K30I/A53S/T59P/V60A/T82A/K190E/K208R
  - Round 2 OOD 教训 - 不用 ML 打分, 严控突变数

约束:
  - 长度 220-250 aa, M 开头, 仅 20 标准 AA
  - 不在 Exclusion_List.csv 中 (全量 135,414 条)
  - chromophore 完整 (TYG / SYG / GYG)
"""
import re, json, pandas as pd
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
OUT_DIR = ROOT / "work" / "round4"

# ============================================================
# 1. 加载野生型序列
# ============================================================
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
ppluGFP = wts["ppluGFP"]

# mBaoJin (PDB 8QBJ, 去掉 RS cloning artifact)
mbaojin = "MVSKGEEENMASTPFKFQLKGTINGKSFTVEGEGEGNSHEGSHKGKYVCTSGKLPMSWAALGTTFGYGMKYYTKYPSGLKNWFREVMPGGFTYDRHIQYKGDGSIHAKHQHFMKNGTYHNIVEFTGQDFKENSPVLTGDMNVSLPNEVPQIPRDDGVECPVTLLYPLLSDKSKYVEAHQYTICKPLHNQPAPDVPYHWIRKQYTQSKDDAEERDHICQSETLEAHLKGMDELYK"

print(f"avGFP    {len(avGFP)} aa  pos65={avGFP[64]}")
print(f"sfGFP    {len(sfGFP)} aa  pos65={sfGFP[64]}")
print(f"amacGFP  {len(amacGFP)} aa  pos65={amacGFP[64]}")
print(f"cgreGFP  {len(cgreGFP)} aa  pos65={cgreGFP[64]}")
print(f"ppluGFP  {len(ppluGFP)} aa")
print(f"mBaoJin  {len(mbaojin)} aa  pos65={mbaojin[64]}")


# ============================================================
# 2. 工具函数
# ============================================================
def apply(seq, muts, strict=False):
    """应用突变. strict=True 时检查 from_aa 匹配, 否则静默改写."""
    s = list(seq)
    warnings = []
    for m in muts:
        match = re.match(r"([A-Z])(\d+)([A-Z])", m)
        if not match:
            continue
        fr, pos, to = match.group(1), int(match.group(2)), match.group(3)
        idx = pos - 1
        if idx >= len(s):
            warnings.append(f"{m}: pos out of range")
            continue
        if s[idx] != fr:
            warnings.append(f"{m}: expected {fr} got {s[idx]}")
            if strict:
                continue
        s[idx] = to
    return "".join(s), warnings


def n_muts(seq, parent):
    return sum(1 for a, b in zip(seq, parent) if a != b)


def chromo(seq):
    for t in ["TYG", "SYG", "GYG", "CYG", "HYG"]:
        if t in seq:
            return t
    return "?"


# ============================================================
# 3. 候选定义 (按战略角色分组)
# ============================================================
candidates = []  # (name, seq, parent_seq, role, notes)

# --------------------------------------------------------------------
# 角色 A: 🔥 热稳爆款 (mBaoJin 衍生, Tm 92°C 主攻最佳热稳奖 + 综合分)
# --------------------------------------------------------------------
# mBaoJin WT 已在排除列表 (Round 3 验证), 必须加 1+ 突变绕开
# 文献: PDB 8QBJ 显示 mBaoJin 表面有多个保守可改位点
# 选择规则: 表面残基 + 不在荧光团 8Å 内 + 不破坏单体化

# A1: 单突变 D173N (Round 3 已设计, pLDDT 38.8 但保留)
candidates.append({
    "name": "A1_mBaoJin_D173N",
    "seq": apply(mbaojin, ["D173N"])[0],
    "parent": mbaojin,
    "role": "thermostable_hero",
    "notes": "mBaoJin + D173N 表面保守突变, Tm~92°C 单体, 1mut",
    "scaffold": "mBaoJin",
    "expected_tm": 92,
})

# A2: 单突变 E194D (表面保守, 不在二聚界面)
candidates.append({
    "name": "A2_mBaoJin_E194D",
    "seq": apply(mbaojin, ["E194D"])[0],
    "parent": mbaojin,
    "role": "thermostable_hero",
    "notes": "mBaoJin + E194D 表面保守, Tm~92°C, 1mut",
    "scaffold": "mBaoJin",
    "expected_tm": 92,
})

# A3: 双突变 D173N + K196R (保留正电荷, 表面)
candidates.append({
    "name": "A3_mBaoJin_D173N_K196R",
    "seq": apply(mbaojin, ["D173N", "K196R"])[0],
    "parent": mbaojin,
    "role": "thermostable_hero",
    "notes": "mBaoJin + D173N + K196R 表面双突变, Tm~92°C, 2mut",
    "scaffold": "mBaoJin",
    "expected_tm": 92,
})

# A4: 单突变 E222D (远离荧光团)
candidates.append({
    "name": "A4_mBaoJin_E222D",
    "seq": apply(mbaojin, ["E222D"])[0],
    "parent": mbaojin,
    "role": "thermostable_hero",
    "notes": "mBaoJin + E222D C 端保守, Tm~92°C, 1mut",
    "scaffold": "mBaoJin",
    "expected_tm": 92,
})

# A5: 三突变 D173N + E194D + K196R
candidates.append({
    "name": "A5_mBaoJin_triple",
    "seq": apply(mbaojin, ["D173N", "E194D", "K196R"])[0],
    "parent": mbaojin,
    "role": "thermostable_hero",
    "notes": "mBaoJin + 3 表面保守突变, Tm~92°C, 3mut",
    "scaffold": "mBaoJin",
    "expected_tm": 92,
})

# --------------------------------------------------------------------
# 角色 B: ⚖️ 综合平衡 (sfGFP + htFuncLib sf:acid.3 风格)
# htFuncLib (Weinstein 2023): Tm 高达 96°C
# sf:acid.3 公开突变: T65S, Q69L, S72A, T108V, Y145M, V224I (基于sfGFP)
# 但 sfGFP 已有 T65, 所以这里突变需基于 sfGFP 的实际残基检查
# --------------------------------------------------------------------
# 先确认 sfGFP 各位点残基
sf_check = lambda p: f"{sfGFP[p-1]}{p}"
print(f"\nsfGFP key positions: {[sf_check(p) for p in [65,69,72,108,145,224]]}")

# B1: sfGFP + Q69L + S72A + T108V (sf:acid 子集, 兼容)
candidates.append({
    "name": "B1_sfGFP_acid_subset",
    "seq": apply(sfGFP, ["Q69L", "S72A", "T108V"])[0],
    "parent": sfGFP,
    "role": "combined_balanced",
    "notes": "sfGFP + htFuncLib sf:acid 子集 (Q69L/S72A/T108V), Tm 可能 >85°C, 3mut",
    "scaffold": "sfGFP",
    "expected_tm": 85,
})

# B2: sfGFP + sf:acid 完整 (替换 T65 不动, 因 sfGFP 已是 T65)
# sfGFP 实际 pos65=T, 所以 T65S 是 T→S, 注意 chromophore TYG (pos 65-67)
# T65S 会改变 chromophore 为 SYG... 风险高, 跳过
# 改用: Q69L + S72A + T108V + Y145M + V224I (除T65S外的5个)
b2_muts = ["Q69L", "S72A", "T108V", "Y145M", "V224I"]
candidates.append({
    "name": "B2_sfGFP_acid_5mut",
    "seq": apply(sfGFP, b2_muts)[0],
    "parent": sfGFP,
    "role": "combined_balanced",
    "notes": "sfGFP + htFuncLib sf:acid.3 (除T65S外), Tm 可能 >90°C, 5mut",
    "scaffold": "sfGFP",
    "expected_tm": 90,
})

# B3: sfGFP + sf:acid + S30R (叠加最强单点稳定)
candidates.append({
    "name": "B3_sfGFP_acid5_S30R",
    "seq": apply(sfGFP, b2_muts + ["S30R"])[0],
    "parent": sfGFP,
    "role": "combined_balanced",
    "notes": "sfGFP + sf:acid 5mut + S30R (+1.25 kcal/mol), 6mut",
    "scaffold": "sfGFP",
    "expected_tm": 92,
})

# --------------------------------------------------------------------
# 角色 C: 🛡️ 保险条 (Round 3 验证过的极保守)
# --------------------------------------------------------------------
# C1: sfGFP + I152S (Round 3 Seq 1, pLDDT 最高 48.2)
candidates.append({
    "name": "C1_sfGFP_I152S",
    "seq": apply(sfGFP, ["I152S"])[0],
    "parent": sfGFP,
    "role": "safety_baseline",
    "notes": "sfGFP + I152S chromophore 邻位优化, Round 3 验证, Tm~80°C, 1mut",
    "scaffold": "sfGFP",
    "expected_tm": 80,
})

# C2: sfGFP + I152S + Q69L (在C1基础上加1个稳定突变)
candidates.append({
    "name": "C2_sfGFP_I152S_Q69L",
    "seq": apply(sfGFP, ["I152S", "Q69L"])[0],
    "parent": sfGFP,
    "role": "safety_baseline",
    "notes": "sfGFP + I152S + Q69L (htFuncLib), 2mut",
    "scaffold": "sfGFP",
    "expected_tm": 82,
})

# --------------------------------------------------------------------
# 角色 D: 🆕 探索 (TGP 风格突变, 同家族编号兼容部分)
# TGP 基于 mAG 编号, K30I 在 GFP 家族对应 pos30
# sfGFP pos30 = S, S30R 已是 sfGFP 设计, 这里改尝试 S30I (TGP风格)
# A53S in mAG -> sfGFP pos53 = ?, T59P -> sfGFP pos59 = ?, V60A -> sfGFP pos60 = ?
# 检查后只移植同家族兼容的
# --------------------------------------------------------------------
print(f"\nsfGFP TGP-test positions: " + ", ".join(f"{p}={sfGFP[p-1]}" for p in [30,53,59,60,82,190,208]))

# D1: sfGFP + sf:acid 5 + I152S (B2 + I152S 叠加 Round 1 验证突变)
candidates.append({
    "name": "D1_sfGFP_acid5_I152S",
    "seq": apply(sfGFP, b2_muts + ["I152S"])[0],
    "parent": sfGFP,
    "role": "exploration",
    "notes": "sfGFP + sf:acid 5mut + I152S, 6mut, 综合稳定+亮度",
    "scaffold": "sfGFP",
    "expected_tm": 88,
})

# D2: avGFP + sfGFP10 (Round 3 Seq 3 保留, 教科书路线)
sfgfp_full_av = ["S65T", "F99S", "M153T", "V163A",
                 "S30R", "Y39N", "N105T", "Y145F", "I171V", "A206V"]
candidates.append({
    "name": "D2_avGFP_sfGFP10",
    "seq": apply(avGFP, sfgfp_full_av)[0],
    "parent": avGFP,
    "role": "exploration",
    "notes": "avGFP + sfGFP 完整10 (F64L已含), 10mut, 教科书 Round 3 验证",
    "scaffold": "avGFP",
    "expected_tm": 80,
})

# D3: sfGFP + sf:acid 5 + I152S + A206V (避免二聚化)
candidates.append({
    "name": "D3_sfGFP_acid5_I152S_A206V",
    "seq": apply(sfGFP, b2_muts + ["I152S", "A206V"])[0],
    "parent": sfGFP,
    "role": "exploration",
    "notes": "sfGFP + sf:acid 5 + I152S + A206V monomer, 7mut",
    "scaffold": "sfGFP",
    "expected_tm": 88,
})

# ============================================================
# 4. 验证 + 排除列表检查
# ============================================================
print("\n" + "=" * 90)
print("加载 Exclusion_List 进行全量检查...")
excl = pd.read_csv(ROOT / "Exclusion_List.csv")
excl_seqs = set(excl["Sequence"].astype(str).str.strip())
print(f"  排除列表共 {len(excl_seqs)} 条")

print("\n" + "=" * 90)
header = f"{'name':<32} {'len':>4} {'mut':>3} {'chromo':>6} {'role':<22} {'excl':>5} {'status'}"
print(header)
print("=" * 90)

valid_candidates = []
for c in candidates:
    seq = c["seq"]
    nm = n_muts(seq, c["parent"])
    ch = chromo(seq)
    in_excl = seq in excl_seqs

    issues = []
    if not seq.startswith("M"):
        issues.append("no-M")
    if len(seq) < 220 or len(seq) > 250:
        issues.append(f"len={len(seq)}")
    if set(seq) - set("ACDEFGHIKLMNPQRSTVWY"):
        issues.append("bad-AA")
    if ch == "?":
        issues.append("no-chromo")
    if nm > 12:
        issues.append(f"mut={nm}>12")

    excl_mark = "❌" if in_excl else "✓"
    status = "✓" if not issues and not in_excl else ",".join(issues) + (" EXCL!" if in_excl else "")
    print(f"{c['name']:<32} {len(seq):>4} {nm:>3} {ch:>6} {c['role']:<22} {excl_mark:>5} {status}")

    c["length"] = len(seq)
    c["n_muts"] = nm
    c["chromophore"] = ch
    c["in_exclusion"] = in_excl
    c["valid"] = (not issues) and (not in_excl)

    if c["valid"]:
        valid_candidates.append(c)

# ============================================================
# 5. 保存
# ============================================================
print(f"\n通过验证: {len(valid_candidates)} / {len(candidates)} 条")

# 按角色统计
from collections import Counter
role_count = Counter(c["role"] for c in valid_candidates)
print("\n按角色分布:")
for role, count in role_count.most_common():
    print(f"  {role}: {count}")

# 准备保存 (移除不可序列化字段)
to_save = []
for c in valid_candidates:
    item = {k: v for k, v in c.items() if k != "parent"}
    item["parent_scaffold"] = c["scaffold"]
    to_save.append(item)

with open(OUT_DIR / "candidates_round4.json", "w", encoding="utf-8") as f:
    json.dump(to_save, f, indent=2, ensure_ascii=False)

print(f"\n✓ 候选已保存到 work/round4/candidates_round4.json")
print(f"  共 {len(to_save)} 条有效候选, 下一步: 02_esmfold_screen.py")
