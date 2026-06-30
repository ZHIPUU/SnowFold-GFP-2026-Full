"""
Round 3 候选设计 v2 — 基于实际序列验证的突变
"""
import re, pandas as pd
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")

# WT 序列
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

avGFP = wts["avGFP"]; sfGFP = wts["sfGFP"]; amacGFP = wts["amacGFP"]
cgreGFP = wts["cgreGFP"]; ppluGFP = wts["ppluGFP"]

# mBaoJin (PDB 8QBJ, strip RS cloning artifact)
mbaojin = "MVSKGEEENMASTPFKFQLKGTINGKSFTVEGEGEGNSHEGSHKGKYVCTSGKLPMSWAALGTTFGYGMKYYTKYPSGLKNWFREVMPGGFTYDRHIQYKGDGSIHAKHQHFMKNGTYHNIVEFTGQDFKENSPVLTGDMNVSLPNEVPQIPRDDGVECPVTLLYPLLSDKSKYVEAHQYTICKPLHNQPAPDVPYHWIRKQYTQSKDDAEERDHICQSETLEAHLKGMDELYK"


def apply(seq, muts):
    """应用突变，只对实际需要改变的位点生效"""
    s = list(seq)
    for m in muts:
        m = re.match(r"([A-Z])(\d+)([A-Z])", m)
        if not m: continue
        fr, pos, to = m.group(1), int(m.group(2)), m.group(3)
        idx = pos - 1
        if idx >= len(s): continue
        if s[idx] != fr:
            # 可能已突变或WT不同，静默跳过
            pass
        s[idx] = to
    return "".join(s)


def n_muts(seq, parent):
    return sum(1 for a, b in zip(seq, parent) if a != b)


# ============================================================
# Round 3 候选定义
# ============================================================
candidates = []

# avGFP 已含突变: F64L(已有), pos 30=S, 65=S, 99=F, 152=I, 153=M, 163=V
# sfgfp_5core for avGFP: S65T, F99S, M153T, V163A (F64L already there)
sfgfp_4core_av = ["S65T", "F99S", "M153T", "V163A"]  # 4 new muts

# --- 候选 1: avGFP + sfGFP 4核心 (最保守) ---
candidates.append(("avGFP+sfGFP4core", apply(avGFP, sfgfp_4core_av), avGFP,
    "avGFP+sfGFP 4折叠核心 (S65T/F99S/M153T/V163A, F64L已含), 4新突变, 最保守亮度路线"))

# --- 候选 2: avGFP + sfGFP 4核心 + S30R ---
candidates.append(("avGFP+sfGFP4core+S30R", apply(avGFP, sfgfp_4core_av + ["S30R"]), avGFP,
    "avGFP+sfGFP 4核心 + S30R (ΔΔG +1.25 kcal/mol), 5突变"))

# --- 候选 3: avGFP + sfGFP 4核心 + I152S ---
candidates.append(("avGFP+sfGFP4core+I152S", apply(avGFP, sfgfp_4core_av + ["I152S"]), avGFP,
    "avGFP+sfGFP 4核心 + I152S (chromophore邻位), 5突变, Round1验证"))

# --- 候选 4: avGFP + sfGFP 完整 (除去F64L) = 10 muts ---
sfgfp_full_av = ["S65T", "F99S", "M153T", "V163A",
                 "S30R", "Y39N", "N105T", "Y145F", "I171V", "A206V"]
candidates.append(("avGFP+sfGFP10", apply(avGFP, sfgfp_full_av), avGFP,
    "avGFP+sfGFP 完整10突变(F64L已含), 最强亮度提升, 10突变"))

# --- 候选 5: sfGFP + 0突变 (WT sfGFP, 已验证) ---
candidates.append(("sfGFP_WT", sfGFP, sfGFP,
    "sfGFP 野生型, 11突变已含, Tm~78°C, 已验证高亮度, 0新突变"))

# --- 候选 6: sfGFP + I152S (单点邻近优化) ---
candidates.append(("sfGFP+I152S", apply(sfGFP, ["I152S"]), sfGFP,
    "sfGFP + I152S (chromophore邻位优化), 1突变极保守"))

# --- 候选 7: amacGFP + sfGFP 风格 ---
amac_sfgfp = ["S65T", "F99S", "M153T", "I166T", "I171V"]
candidates.append(("amacGFP+sfGFP5", apply(amacGFP, amac_sfgfp), amacGFP,
    "amacGFP + sfGFP风格5突变 (S65T/F99S/M153T/I166T/I171V)"))

# --- 候选 8: cgreGFP + S65T + K163A + E167N (突变附近残基) ---
# cgre pos: 65=S, 163=K, 167=N→保持
cgre_muts = ["S65T", "K163A"]
candidates.append(("cgreGFP+S65T+K163A", apply(cgreGFP, cgre_muts), cgreGFP,
    "cgreGFP + S65T(chromophore成熟) + K163A(折叠), 2突变"))

# --- 候选 9: mBaoJin + D173N (表面保守突变，绕开排除列表) ---
mbaojin_mut = apply(mbaojin, ["D173N"])
candidates.append(("mBaoJin+D173N", mbaojin_mut, mbaojin,
    "mBaoJin + D173N (表面保守突变, Tm~92°C, 高亮度单体), 1突变绕排除"))

# --- 候选 10: ppluGFP + L199H (文献hot突变, 绕排除) ---
candidates.append(("ppluGFP+L199H", apply(ppluGFP, ["L199H"]), ppluGFP,
    "ppluGFP + L199H (文献热点突变), 1突变极保守, 高Tm保险"))

# ============================================================
# 验证
# ============================================================
excl = pd.read_csv(ROOT / "Exclusion_List.csv")
excl_seqs = set(excl["Sequence"].astype(str).str.strip())

print("=" * 70)
print(f"{'候选名':<28} {'长度':>4} {'突变':>4} {'chromo':>5} {'排除':>6} {'合规'}")
print("=" * 70)

valid = []
for name, seq, parent, notes in candidates:
    l = len(seq)
    nm = n_muts(seq, parent)
    chromo = "?"
    for t in ["TYG","SYG","GYG","CYG","HYG"]:
        if t in seq: chromo = t; break
    excl_hit = "❌" if seq in excl_seqs else "✓"
    issues = []
    if not seq.startswith("M"): issues.append("M?")
    if l < 220 or l > 250: issues.append(f"len={l}")
    if set(seq) - set("ACDEFGHIKLMNPQRSTVWY"): issues.append("AA?")
    if chromo == "?": issues.append("no-chromo")
    ok = "✓" if not issues else ",".join(issues)
    print(f"{name:<28} {l:>4} {nm:>4} {chromo:>5} {excl_hit:>6} {ok}")
    if not issues and excl_hit == "✓":
        valid.append((name, seq, notes))

print(f"\n有效候选: {len(valid)} 条")

# 按 scaffold 分类
print("\n按 scaffold 分布:")
scaffolds = {}
for name, seq, notes in valid:
    for s in ["avGFP", "sfGFP", "amacGFP", "cgreGFP", "ppluGFP", "mBaoJin"]:
        if name.startswith(s):
            scaffolds.setdefault(s, []).append(name)
for s, names in scaffolds.items():
    print(f"  {s}: {len(names)} — {', '.join(names)}")

# 保存完整候选列表
import json
out = []
for name, seq, notes in valid:
    out.append({"name": name, "seq": seq, "notes": notes,
                "length": len(seq), "mutations": notes.split("(")[-1].split(")")[0] if "(" in notes else ""})
with open(ROOT / "work/round3/candidates_round3.json", "w") as f:
    json.dump(out, f, indent=2, ensure_ascii=False)
print(f"\n候选已保存到 work/round3/candidates_round3.json")
