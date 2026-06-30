"""
Round 3 候选设计
策略: 多 scaffold 覆盖 + 控制突变数 ≤10 + 论文知识驱动

设计原则:
1. 每候选突变数控制在 ≤10 (Arcadia 2025: 汉明距离>4 功能急剧下降)
2. 4+ scaffold 覆盖 (avGFP, sfGFP, amacGFP, cgreGFP, mBaoJin)
3. 每候选 1 个核心机制 (不堆砌)
4. 必须 chromophore 完整
"""
import re
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")

# ============================================================
# WT 序列
# ============================================================
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

avGFP = wts["avGFP"]
sfGFP = wts["sfGFP"]
amacGFP = wts["amacGFP"]
cgreGFP = wts["cgreGFP"]
ppluGFP = wts["ppluGFP"]

# mBaoJin (PDB 8QBJ, 去除RS前缀)
mbaojin = "MVSKGEEENMASTPFKFQLKGTINGKSFTVEGEGEGNSHEGSHKGKYVCTSGKLPMSWAALGTTFGYGMKYYTKYPSGLKNWFREVMPGGFTYDRHIQYKGDGSIHAKHQHFMKNGTYHNIVEFTGQDFKENSPVLTGDMNVSLPNEVPQIPRDDGVECPVTLLYPLLSDKSKYVEAHQYTICKPLHNQPAPDVPYHWIRKQYTQSKDDAEERDHICQSETLEAHLKGMDELYK"

print(f"mBaoJin: {len(mbaojin)} aa, chromophore GYG, Tm~92°C")


def apply_muts(parent_seq, mut_list):
    """对 parent 序列应用突变列表, 返回新序列"""
    seq = list(parent_seq)
    for m in mut_list:
        # m format: "S65T" meaning pos 65, from S to T
        match = re.match(r"([A-Z])(\d+)([A-Z])", m)
        if match:
            from_aa, pos, to_aa = match.group(1), int(match.group(2)), match.group(3)
            idx = pos - 1  # 0-based
            if idx < len(seq) and seq[idx] == from_aa:
                seq[idx] = to_aa
            else:
                print(f"  WARNING: {m} — expected {from_aa} at pos {pos}, got {seq[idx] if idx < len(seq) else 'OOB'}")
        else:
            print(f"  WARNING: Can't parse mutation {m}")
    return "".join(seq)


def verify_candidate(name, seq, parent, design_notes):
    """验证候选"""
    issues = []
    if not seq.startswith("M"):
        issues.append("不以M开头")
    if len(seq) < 220 or len(seq) > 250:
        issues.append(f"长度{len(seq)}")
    invalid = {c for c in seq if c not in "ACDEFGHIKLMNPQRSTVWY"}
    if invalid:
        issues.append(f"非标准AA:{invalid}")
    chromo = None
    for tri in ["TYG", "SYG", "GYG", "CYG", "HYG"]:
        if tri in seq:
            chromo = tri
            break
    if not chromo:
        issues.append("无chromophore")

    status = "✓" if not issues else "✗ " + " | ".join(issues)
    n_mut = sum(1 for a, b in zip(seq, parent) if a != b)
    print(f"  {name}: {len(seq)}aa, {n_mut} muts, chromophore={chromo}  {status}")
    print(f"    {design_notes}")
    return len(issues) == 0


# ============================================================
# 候选设计
# ============================================================
candidates = {}

# --- 候选 A: avGFP + sfGFP 5 折叠核心 (最安全) ---
# 5 muts: F64L/S65T/F99S/M153T/V163A
name_a = "avGFP+sfGFP5core"
seq_a = apply_muts(avGFP, ["F64L", "S65T", "F99S", "M153T", "V163A"])
candidates[name_a] = (seq_a, avGFP, "avGFP + sfGFP 折叠核心5突变 (F64L/S65T/F99S/M153T/V163A), 已验证的保守路线")

# --- 候选 B: avGFP + sfGFP 5核心 + S30R (亮度增强) ---
# 6 muts: core5 + S30R
name_b = "avGFP+sfGFP5core+S30R"
seq_b = apply_muts(avGFP, ["F64L", "S65T", "F99S", "M153T", "V163A", "S30R"])
candidates[name_b] = (seq_b, avGFP, "avGFP + sfGFP 5核心 + S30R (+1.25 kcal/mol), 6突变稳健增强")

# --- 候选 C: avGFP + sfGFP 5核心 + I152S (Round1主力) ---
name_c = "avGFP+sfGFP5core+I152S"
seq_c = apply_muts(avGFP, ["F64L", "S65T", "F99S", "M153T", "V163A", "I152S"])
candidates[name_c] = (seq_c, avGFP, "avGFP + sfGFP 5核心 + I152S (chromophore邻位), Round1验证有效 → 6突变")

# --- 候选 D: avGFP + sfGFP 完整 11 突变 (标准sfGFP化) ---
name_d = "avGFP+sfGFP11"
sfgfp_muts = ["F64L", "S65T", "F99S", "M153T", "V163A",
              "S30R", "Y39N", "N105T", "Y145F", "I171V", "A206V"]
seq_d = apply_muts(avGFP, sfgfp_muts)
candidates[name_d] = (seq_d, avGFP, "avGFP + sfGFP 完整11突变, 最强已验证的亮度提升路线")

# --- 候选 E: avGFP + sfGFP 5核心 + TGP 稳定核心 (4个) = 9 muts ---
# TGP 4 稳定: A53S/T59P/V60A/T82A
name_e = "avGFP+sfGFP5+TGP4stable"
tgp_stable_avGFP = ["A53S", "T59P", "V60A", "T82A"]
seq_e = apply_muts(avGFP, ["F64L", "S65T", "F99S", "M153T", "V163A"] + tgp_stable_avGFP)
candidates[name_e] = (seq_e, avGFP, "avGFP + sfGFP 5核心 + TGP 4稳定核心 (A53S/T59P/V60A/T82A), 9突变")

# --- 候选 F: mBaoJin 原始 (超高Tm保险条) ---
name_f = "mBaoJin_WT"
candidates[name_f] = (mbaojin, mbaojin, "mBaoJin 原始 (单体StayGold), Tm~92°C, 高亮度, 快速成熟 — 最强保险条")

# --- 候选 G: cgreGFP + S65T (最小改动) ---
name_g = "cgreGFP+S65T"
seq_g = apply_muts(cgreGFP, ["S65T"])
candidates[name_g] = (seq_g, cgreGFP, "cgreGFP + S65T (chromophore成熟加速), 1突变极保守, baseline 4.50高")

# --- 候选 H: sfGFP + TGP 稳定核心 (4 muts) ---
name_h = "sfGFP+TGP4stable"
seq_h = apply_muts(sfGFP, ["A53S", "T59P", "V60A", "T82A"])
candidates[name_h] = (seq_h, sfGFP, "sfGFP + TGP 4稳定核心, 4突变极保守, sfGFP已验证骨架")

# --- 候选 I: amacGFP + sfGFP 风格 (6 muts) ---
name_i = "amacGFP+sfGFPstyle6"
amac_muts = ["S65T", "V99S", "M153T", "I166T", "I171V", "A206V"]
seq_i = apply_muts(amacGFP, amac_muts)
candidates[name_i] = (seq_i, amacGFP, "amacGFP + sfGFP风格6突变 (S65T/V99S/M153T/I166T/I171V/A206V)")

# ============================================================
# 验证所有候选
# ============================================================
print("\n" + "=" * 60)
print("候选验证")
print("=" * 60)

valid_candidates = {}
for name, (seq, parent, notes) in candidates.items():
    ok = verify_candidate(name, seq, parent, notes)
    if ok:
        valid_candidates[name] = (seq, notes)

# ============================================================
# 排除列表检查
# ============================================================
print("\n" + "=" * 60)
print("排除列表检查")
print("=" * 60)
import pandas as pd
excl = pd.read_csv(ROOT / "Exclusion_List.csv")
excl_seqs = set(excl["Sequence"].astype(str).str.strip())

for name, (seq, notes) in valid_candidates.items():
    in_excl = seq in excl_seqs
    status = "❌ 排除列表!" if in_excl else "✓"
    print(f"  {name}: {status}")

# ============================================================
# 总结
# ============================================================
print("\n" + "=" * 60)
print("Round 3 候选汇总")
print("=" * 60)

# 按 scaffold 分组
scaffolds = {"avGFP": [], "sfGFP": [], "amacGFP": [], "cgreGFP": [], "mBaoJin": []}
for name, (seq, notes) in valid_candidates.items():
    if "avGFP" in name:
        scaffolds["avGFP"].append(name)
    elif "sfGFP" in name:
        scaffolds["sfGFP"].append(name)
    elif "amacGFP" in name:
        scaffolds["amacGFP"].append(name)
    elif "cgreGFP" in name:
        scaffolds["cgreGFP"].append(name)
    elif "mBaoJin" in name:
        scaffolds["mBaoJin"].append(name)

for scaff, names in scaffolds.items():
    if names:
        print(f"  {scaff}: {len(names)} candidates — {', '.join(names)}")

print(f"\n  总计 {len(valid_candidates)} 条有效候选")
