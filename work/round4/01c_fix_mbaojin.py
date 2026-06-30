"""
Round 4 Step 1C: 修正 mBaoJin 候选 (基于真实残基)
==================================
之前错误: 假设 mBaoJin pos 173=D, 194=E, 222=E
实际:           pos 173=K, 194=V, 222=L

正确策略:
  - 表面残基突变 (远离 chromophore pos 63-65 GYG)
  - 不破坏单体化关键位点 (Q140P, H141Q, C165Y, N171Y from mBaoJin paper)
  - 选择保守替换 (K->R, E->D, V->I 等)
"""
import re, json, pandas as pd
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")
OUT = ROOT / "work" / "round4"

mbaojin = "MVSKGEEENMASTPFKFQLKGTINGKSFTVEGEGEGNSHEGSHKGKYVCTSGKLPMSWAALGTTFGYGMKYYTKYPSGLKNWFREVMPGGFTYDRHIQYKGDGSIHAKHQHFMKNGTYHNIVEFTGQDFKENSPVLTGDMNVSLPNEVPQIPRDDGVECPVTLLYPLLSDKSKYVEAHQYTICKPLHNQPAPDVPYHWIRKQYTQSKDDAEERDHICQSETLEAHLKGMDELYK"

print(f"mBaoJin length: {len(mbaojin)}")
print(f"Chromophore region (60-68): {mbaojin[59:68]}")
# 找出 GYG (chromophore)
for i in range(len(mbaojin)-2):
    if mbaojin[i:i+3] == "GYG":
        print(f"GYG chromophore at pos {i+1}-{i+3}")

# 显示完整序列, 每行 60 个字符
print("\n完整 mBaoJin 序列:")
for i in range(0, len(mbaojin), 60):
    nums = "".join(f"{(i+j+1)%10}" for j in range(min(60, len(mbaojin)-i)))
    print(f"  {i+1:>4}: {nums}")
    print(f"        {mbaojin[i:i+60]}")

# ============================================================
# 真实表面残基 (远离 chromophore, 表面暴露)
# ============================================================
# 选择策略:
#   1. C 端柔性区 (pos 220-234, 远离 GYG)
#   2. 表面 loop 区域 (pos 100-110, 140-160 之间且不在单体化关键位点)
#   3. 保守替换 K↔R, E↔D, V↔I, S↔T, Q↔N
# ============================================================

# 真实可用的突变 (基于序列扫描):
print("\nmBaoJin 候选突变扫描:")
candidate_muts = []
for p in range(150, 234):  # 远离 chromophore
    aa = mbaojin[p-1]
    # 保守替换映射
    cons = {"K": "R", "R": "K", "E": "D", "D": "E", "S": "T", "T": "S",
            "V": "I", "I": "V", "Q": "N", "N": "Q", "L": "M", "M": "L"}
    if aa in cons:
        target = cons[aa]
        m = f"{aa}{p}{target}"
        candidate_muts.append((p, m))

# 选择不在 mBaoJin 论文关键单体化位点 (PDB编号S55T/H77R/E80G/Q140P/H141Q/C165Y/N171Y/T201A)
# 这些位点 mBaoJin 已突变, 不要回退
forbidden_pos = set()  # 已是 mBaoJin 突变后的残基, 不要动

# 选 10-15 个不同位置的保守突变
selected_singles = []
for p, m in candidate_muts:
    if p < 165 or p > 220:  # 避开 chromophore 周围 + 包含 C 端
        if p in [173, 175, 180, 196, 200, 205, 210, 215, 220, 225, 230]:
            selected_singles.append(m)

print(f"\n选定的单点突变 (基于真实序列):")
for m in selected_singles:
    print(f"  {m}")

# ============================================================
# 应用函数
# ============================================================
def apply_safe(seq, muts, name):
    """应用突变并验证"""
    s = list(seq)
    warnings = []
    for m in muts:
        match = re.match(r"([A-Z])(\d+)([A-Z])", m)
        if not match: continue
        fr, pos, to = match.group(1), int(match.group(2)), match.group(3)
        idx = pos - 1
        if idx >= len(s):
            warnings.append(f"{m}: out of range")
            continue
        if s[idx] != fr:
            warnings.append(f"{m}: WT was {s[idx]}, NOT {fr}")
            continue  # 严格模式: 不应用
        s[idx] = to
    if warnings:
        print(f"  ⚠ {name}: {warnings}")
    return "".join(s)


# ============================================================
# 候选 (基于真实 mBaoJin 残基)
# ============================================================
candidates = []

# 真实可用突变集
real_mbaojin_muts = [
    ("M1_mBaoJin_K173R", ["K173R"], "K173R 保守正电荷, 表面 1mut"),
    ("M2_mBaoJin_V194I", ["V194I"], "V194I 疏水保守, 1mut"),
    ("M3_mBaoJin_L222M", ["L222M"], "L222M 疏水保守, C端 1mut"),
    ("M4_mBaoJin_Y196F", ["Y196F"], "Y196F 芳香保守, 1mut"),
    ("M5_mBaoJin_D230E", ["D230E"], "D230E 同电荷, C端 1mut"),
    ("M6_mBaoJin_K173R_V194I", ["K173R", "V194I"], "K173R + V194I 双保守, 2mut"),
    ("M7_mBaoJin_K173R_Y196F", ["K173R", "Y196F"], "K173R + Y196F 双保守, 2mut"),
    ("M8_mBaoJin_combo3", ["K173R", "V194I", "L222M"], "K173R + V194I + L222M 三保守, 3mut"),
    ("M9_mBaoJin_E180D", ["E180D"], "E180D 同电荷 (pos180 实际值待检)"),
    ("M10_mBaoJin_K200R", ["K200R"], "pos200 检查"),
]

for name, muts, notes in real_mbaojin_muts:
    seq = apply_safe(mbaojin, muts, name)
    candidates.append({
        "name": name,
        "seq": seq,
        "parent": mbaojin,
        "role": "thermostable_hero",
        "notes": notes,
        "scaffold": "mBaoJin",
        "expected_tm": 92,
    })


# ============================================================
# 验证
# ============================================================
def n_muts(seq, parent):
    return sum(1 for a, b in zip(seq, parent) if a != b)

def chromo(seq):
    for t in ["TYG", "SYG", "GYG", "CYG", "HYG"]:
        if t in seq: return t
    return "?"

excl_seqs = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())

print("\n" + "=" * 90)
print(f"{'name':<32} {'len':>4} {'mut':>3} {'chromo':>6} {'excl':>5} {'status'}")
print("=" * 90)

valid_new = []
for c in candidates:
    seq = c["seq"]
    nm = n_muts(seq, c["parent"])
    ch = chromo(seq)
    in_excl = seq in excl_seqs

    issues = []
    if not seq.startswith("M"): issues.append("no-M")
    if len(seq) < 220 or len(seq) > 250: issues.append(f"len={len(seq)}")
    if ch == "?": issues.append("no-chromo")
    if nm == 0: issues.append("zero-mut")  # 突变没应用

    excl_mark = "❌" if in_excl else "✓"
    status = "✓" if not issues and not in_excl else ",".join(issues) + (" EXCL!" if in_excl else "")
    print(f"{c['name']:<32} {len(seq):>4} {nm:>3} {ch:>6} {excl_mark:>5} {status}")

    c["length"] = len(seq)
    c["n_muts"] = nm
    c["chromophore"] = ch
    c["in_exclusion"] = in_excl
    c["valid"] = (not issues) and (not in_excl)

    if c["valid"]:
        valid_new.append(c)

print(f"\n通过: {len(valid_new)} / {len(candidates)}")

# 合并到 extended 候选池
with open(OUT / "candidates_round4_extended.json", encoding="utf-8") as f:
    existing = json.load(f)

# 删除旧的错误 A1-A11 mBaoJin 候选 (前缀 A 开头)
filtered = [c for c in existing if not (c["name"].startswith("A") and c["scaffold"] == "mBaoJin")]
print(f"\n移除错误 A1-A11 mBaoJin 候选: {len(existing) - len(filtered)} 条")

# 加入新的 M1-M10
all_seqs = set(c["seq"] for c in filtered)
added = 0
for c in valid_new:
    if c["seq"] in all_seqs: continue
    c_save = {k: v for k, v in c.items() if k != "parent"}
    c_save["parent_scaffold"] = c["scaffold"]
    filtered.append(c_save)
    all_seqs.add(c["seq"])
    added += 1

print(f"新增正确 M1-M10: {added} 条")
print(f"最终候选池: {len(filtered)} 条")

with open(OUT / "candidates_round4_extended.json", "w", encoding="utf-8") as f:
    json.dump(filtered, f, indent=2, ensure_ascii=False)

from collections import Counter
print("\n按角色分布:")
for role, count in Counter(c["role"] for c in filtered).most_common():
    print(f"  {role}: {count}")

print("\n下一步: 重跑 02b_esmfold_screen_extended.py")
