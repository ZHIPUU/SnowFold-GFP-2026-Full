"""
Round 4 Step 1B: 扩展候选池
==================================
基于 01_seed_candidates 的发现:
  - sfGFP 已含 S30R/Y39N/N105T/Y145F/I171V/A206V (所以这些突变在 sfGFP 上是 silent)
  - 需要更多 mBaoJin 变体
  - 需要 avGFP 起点的"无歧义"突变路线 (避免 sfGFP 上 silent mutation 导致重复)
"""
import re, json, pandas as pd
from pathlib import Path

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
mbaojin = "MVSKGEEENMASTPFKFQLKGTINGKSFTVEGEGEGNSHEGSHKGKYVCTSGKLPMSWAALGTTFGYGMKYYTKYPSGLKNWFREVMPGGFTYDRHIQYKGDGSIHAKHQHFMKNGTYHNIVEFTGQDFKENSPVLTGDMNVSLPNEVPQIPRDDGVECPVTLLYPLLSDKSKYVEAHQYTICKPLHNQPAPDVPYHWIRKQYTQSKDDAEERDHICQSETLEAHLKGMDELYK"

def apply(seq, muts):
    s = list(seq)
    for m in muts:
        match = re.match(r"([A-Z])(\d+)([A-Z])", m)
        if not match: continue
        fr, pos, to = match.group(1), int(match.group(2)), match.group(3)
        idx = pos - 1
        if idx >= len(s): continue
        s[idx] = to
    return "".join(s)

def n_muts(seq, parent):
    return sum(1 for a, b in zip(seq, parent) if a != b)

def chromo(seq):
    for t in ["TYG", "SYG", "GYG", "CYG", "HYG"]:
        if t in seq: return t
    return "?"

# 检查 mBaoJin 表面残基 (远离 chromophore pos 60-70)
print("mBaoJin 残基扫描 (后半段表面候选):")
for p in [173, 180, 185, 190, 194, 196, 200, 205, 210, 215, 220, 222, 225, 230]:
    if p <= len(mbaojin):
        print(f"  pos {p}: {mbaojin[p-1]}")

# ============================================================
# 候选定义 (修复重复 + 扩展 mBaoJin)
# ============================================================
candidates = []

# --------------------------------------------------------------------
# 角色 A: 🔥 mBaoJin 系列 (扩展)
# --------------------------------------------------------------------
# 已有 (Round 4-v1):
#   A1: D173N (1mut)
#   A2: E194D (1mut)  
#   A3: D173N + K196R (2mut)
#   A4: E222D (1mut)
#   A5: D173N + E194D + K196R (3mut)
# 新增:
#   A6: D173E (D->E, 保守同电荷)
#   A7: D173N + Q230H (双表面)
#   A8: E194Q + K196R (中性化)
#   A9: D173N + E222D (双 D->N/D)
#   A10: E222Q + Q230H (C 端组合)
new_mbaojin_muts = [
    ("A6_mBaoJin_D173E", ["D173E"], "D173E (保守同电荷), Tm~92, 1mut"),
    ("A7_mBaoJin_D173N_Q230H", ["D173N", "Q230H"], "D173N + Q230H 双表面保守, Tm~92, 2mut"),
    ("A8_mBaoJin_E194Q_K196R", ["E194Q", "K196R"], "E194Q + K196R 中性化, Tm~92, 2mut"),
    ("A9_mBaoJin_D173N_E222D", ["D173N", "E222D"], "D173N + E222D 双 D 突变, Tm~92, 2mut"),
    ("A10_mBaoJin_C_combo", ["E222Q", "Q230H"], "E222Q + Q230H C端组合, Tm~92, 2mut"),
    ("A11_mBaoJin_D173N_E194D", ["D173N", "E194D"], "D173N + E194D 双表面, Tm~92, 2mut"),
]

for name, muts, notes in new_mbaojin_muts:
    seq = apply(mbaojin, muts)
    candidates.append({
        "name": name,
        "seq": seq,
        "parent": mbaojin,
        "role": "thermostable_hero",
        "notes": notes,
        "scaffold": "mBaoJin",
        "expected_tm": 92,
    })

# --------------------------------------------------------------------
# 角色 B: 综合平衡 (基于 avGFP 起点, 避免 sfGFP silent)
# --------------------------------------------------------------------
# avGFP -> sfGFP 11 muts, 加 htFuncLib sf:acid 突变
# avGFP 位置: 30=S, 39=Y, 65=S, 105=N, 145=Y, 171=I, 206=A (野生)
# sfGFP 后的 sf:acid 还可加: Q69L (Q->L), S72A (S->A), T108V (T->V), Y145M
# 但要从 avGFP 起点, 先变成 sfGFP, 再加 sf:acid

# B4: avGFP + sfGFP10 + Q69L (htFuncLib风格的"超 sfGFP")
b4_muts = ["S65T", "F99S", "M153T", "V163A",  # sfGFP 4 core
           "S30R", "Y39N", "N105T", "Y145F", "I171V", "A206V",  # sfGFP 6 表面
           "Q69L"]  # htFuncLib sf:acid 加成
candidates.append({
    "name": "B4_avGFP_sfGFP10_Q69L",
    "seq": apply(avGFP, b4_muts),
    "parent": avGFP,
    "role": "combined_balanced",
    "notes": "avGFP + sfGFP 10 + Q69L htFuncLib, 11mut, 'super-sfGFP'",
    "scaffold": "avGFP",
    "expected_tm": 88,
})

# B5: avGFP + sfGFP10 + Q69L + S72A
b5_muts = b4_muts + ["S72A"]
candidates.append({
    "name": "B5_avGFP_sfGFP10_acid2",
    "seq": apply(avGFP, b5_muts),
    "parent": avGFP,
    "role": "combined_balanced",
    "notes": "avGFP + sfGFP 10 + Q69L + S72A htFuncLib, 12mut",
    "scaffold": "avGFP",
    "expected_tm": 90,
})

# B6: avGFP + sfGFP 4 core + Q69L (轻量 htFuncLib)
b6_muts = ["S65T", "F99S", "M153T", "V163A", "Q69L"]
candidates.append({
    "name": "B6_avGFP_sf4_Q69L",
    "seq": apply(avGFP, b6_muts),
    "parent": avGFP,
    "role": "combined_balanced",
    "notes": "avGFP + sfGFP 4 core + Q69L, 5mut 轻量化路线",
    "scaffold": "avGFP",
    "expected_tm": 80,
})

# --------------------------------------------------------------------
# 角色 C: 保险条 (avGFP 不同变体)
# --------------------------------------------------------------------
# C3: avGFP + sfGFP 4 core + I152S (Round 1 Seq 6 类型, 综合分 9.50 估)
c3_muts = ["S65T", "F99S", "M153T", "V163A", "I152S"]
candidates.append({
    "name": "C3_avGFP_sf4_I152S",
    "seq": apply(avGFP, c3_muts),
    "parent": avGFP,
    "role": "safety_baseline",
    "notes": "avGFP + sfGFP 4 core + I152S, Round 1 综合分 9.50 估计版, 5mut",
    "scaffold": "avGFP",
    "expected_tm": 80,
})

# --------------------------------------------------------------------
# 角色 D: 探索 (TGP 风格 + 多家族)
# --------------------------------------------------------------------
# D4: sfGFP + I152S + K166V (cgreGFP-style hotspot)
candidates.append({
    "name": "D4_sfGFP_I152S_K166V",
    "seq": apply(sfGFP, ["I152S", "K166V"]),
    "parent": sfGFP,
    "role": "exploration",
    "notes": "sfGFP + I152S + K166V cgre-style hotspot, 2mut",
    "scaffold": "sfGFP",
    "expected_tm": 80,
})

# D5: avGFP + sfGFP 4 core + Q69L + I152S (htFuncLib + Round 1)
d5_muts = ["S65T", "F99S", "M153T", "V163A", "Q69L", "I152S"]
candidates.append({
    "name": "D5_avGFP_sf4_acid_I152S",
    "seq": apply(avGFP, d5_muts),
    "parent": avGFP,
    "role": "exploration",
    "notes": "avGFP + sf4 + Q69L + I152S, 6mut 多机制叠加",
    "scaffold": "avGFP",
    "expected_tm": 85,
})

# ============================================================
# 验证 + 排除列表检查
# ============================================================
print("\n" + "=" * 95)
print("加载 Exclusion_List 全量检查...")
excl_seqs = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())
print(f"  排除列表: {len(excl_seqs)} 条")

print("\n" + "=" * 95)
print(f"{'name':<32} {'len':>4} {'mut':>3} {'chromo':>6} {'role':<22} {'excl':>5} {'status'}")
print("=" * 95)

valid_new = []
for c in candidates:
    seq = c["seq"]
    nm = n_muts(seq, c["parent"])
    ch = chromo(seq)
    in_excl = seq in excl_seqs

    issues = []
    if not seq.startswith("M"): issues.append("no-M")
    if len(seq) < 220 or len(seq) > 250: issues.append(f"len={len(seq)}")
    if set(seq) - set("ACDEFGHIKLMNPQRSTVWY"): issues.append("bad-AA")
    if ch == "?": issues.append("no-chromo")

    excl_mark = "❌" if in_excl else "✓"
    status = "✓" if not issues and not in_excl else ",".join(issues) + (" EXCL!" if in_excl else "")
    print(f"{c['name']:<32} {len(seq):>4} {nm:>3} {ch:>6} {c['role']:<22} {excl_mark:>5} {status}")

    c["length"] = len(seq)
    c["n_muts"] = nm
    c["chromophore"] = ch
    c["in_exclusion"] = in_excl
    c["valid"] = (not issues) and (not in_excl)

    if c["valid"]:
        valid_new.append(c)

print(f"\n新增有效候选: {len(valid_new)}")

# ============================================================
# 合并旧候选, 去重
# ============================================================
with open(OUT / "candidates_round4.json", encoding="utf-8") as f:
    old_candidates = json.load(f)

all_seqs = set(c["seq"] for c in old_candidates)
final_candidates = list(old_candidates)
n_added = 0
for c in valid_new:
    if c["seq"] not in all_seqs:
        c_save = {k: v for k, v in c.items() if k != "parent"}
        c_save["parent_scaffold"] = c["scaffold"]
        final_candidates.append(c_save)
        all_seqs.add(c["seq"])
        n_added += 1

print(f"\n合并去重: 旧 {len(old_candidates)} + 新增 {n_added} = 总 {len(final_candidates)}")

# 保存扩展版
with open(OUT / "candidates_round4_extended.json", "w", encoding="utf-8") as f:
    json.dump(final_candidates, f, indent=2, ensure_ascii=False)

print(f"\n✓ 扩展候选池保存到 candidates_round4_extended.json")
print(f"  下一步: 重跑 02_esmfold_screen 评估新候选")

# 按角色统计
from collections import Counter
role_count = Counter(c["role"] for c in final_candidates)
print(f"\n按角色分布:")
for role, count in role_count.most_common():
    print(f"  {role}: {count}")
