"""
Round 4 Step B: 扩展高分预期候选
==================================
基于得分预测发现 (07_score_estimate.py):
  - X4 (avGFP+sfGFP10+I152S+Q69L) 预测综合分 1.36 (最高)
  - sfGFP+Q69L 路线 (htFuncLib Tm 96°C) 未被充分挖掘
  - mBaoJin pLDDT 偏低导致综合分 0.67 偏低

新候选策略 (扩展 X4 公式):
  1. avGFP/sfGFP + sfGFP10 + I152S + Q69L + S72A + T108V (sf:acid 全套)
  2. sfGFP + I152S + Q69L + S72A (轻量 htFuncLib)
  3. sfGFP + Q69L + S72A + T108V (纯 htFuncLib)
  4. avGFP + sfGFP10 + I152S + Q69L + K166V (新组合)
  5. mBaoJin 替代: avGFP+sfGFP+Q69L+I152S 双高分备选
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

def apply(seq, muts):
    s = list(seq)
    for m in muts:
        match = re.match(r"([A-Z])(\d+)([A-Z])", m)
        if not match: continue
        fr, pos, to = match.group(1), int(match.group(2)), match.group(3)
        idx = pos - 1
        if idx >= len(s): continue
        if s[idx] != fr:
            continue
        s[idx] = to
    return "".join(s)

def n_muts(seq, parent):
    return sum(1 for a, b in zip(seq, parent) if a != b)

def chromo(seq):
    for t in ["TYG", "SYG", "GYG"]:
        if t in seq: return t
    return "?"

# ============================================================
# 高分预期新候选
# ============================================================
candidates = []

# 检查 sfGFP/avGFP 关键位点
print(f"sfGFP positions: 69={sfGFP[68]}, 72={sfGFP[71]}, 108={sfGFP[107]}, 145={sfGFP[144]}, 152={sfGFP[151]}, 166={sfGFP[165]}, 224={sfGFP[223]}")
print(f"avGFP positions: 69={avGFP[68]}, 72={avGFP[71]}, 108={avGFP[107]}, 145={avGFP[144]}, 152={avGFP[151]}, 166={avGFP[165]}, 224={avGFP[223]}")

# sfGFP 基础完整突变 (avGFP -> sfGFP 10 mut, F64L已含)
sfgfp_full_av = ["S65T", "F99S", "M153T", "V163A",
                 "S30R", "Y39N", "N105T", "Y145F", "I171V", "A206V"]

# --------------------------------------------------------------------
# H1-H5: avGFP 系列升级 (X4 路线扩展)
# --------------------------------------------------------------------
# H1: avGFP + sfGFP10 + I152S + Q69L + S72A (X4 + S72A)
candidates.append({
    "name": "H1_avGFP_sfGFP_acid3_I152S",
    "seq": apply(avGFP, sfgfp_full_av + ["I152S", "Q69L", "S72A"]),
    "parent": avGFP, "scaffold": "avGFP",
    "role": "combined_hero", "expected_tm": 92,
    "notes": "avGFP + sfGFP10 + I152S + Q69L + S72A 强力组合, 13mut"
})

# H2: avGFP + sfGFP10 + I152S + Q69L + T108V (X4 + T108V)
candidates.append({
    "name": "H2_avGFP_sfGFP_Q69L_T108V",
    "seq": apply(avGFP, sfgfp_full_av + ["I152S", "Q69L", "T108V"]),
    "parent": avGFP, "scaffold": "avGFP",
    "role": "combined_hero", "expected_tm": 90,
    "notes": "avGFP + sfGFP10 + I152S + Q69L + T108V, 13mut"
})

# H3: avGFP + sfGFP10 + I152S + Q69L + S72A + T108V (X4 + 完整htFuncLib)
candidates.append({
    "name": "H3_avGFP_sfGFP_acid_full",
    "seq": apply(avGFP, sfgfp_full_av + ["I152S", "Q69L", "S72A", "T108V"]),
    "parent": avGFP, "scaffold": "avGFP",
    "role": "combined_hero", "expected_tm": 94,
    "notes": "avGFP + sfGFP10 + I152S + 完整htFuncLib(Q69L/S72A/T108V), 14mut"
})

# H4: avGFP + sfGFP10 + I152S + Q69L + K166V (X4 + K166V cgre hotspot)
candidates.append({
    "name": "H4_avGFP_sfGFP_Q69L_K166V",
    "seq": apply(avGFP, sfgfp_full_av + ["I152S", "Q69L", "K166V"]),
    "parent": avGFP, "scaffold": "avGFP",
    "role": "combined_hero", "expected_tm": 90,
    "notes": "avGFP + sfGFP10 + I152S + Q69L + K166V cgre-hotspot, 13mut"
})

# H5: avGFP + sfGFP10 (Round 3 D2) + Q69L (sfGFP 上 silent, avGFP 上有效)
candidates.append({
    "name": "H5_avGFP_sfGFP10_Q69L",
    "seq": apply(avGFP, sfgfp_full_av + ["Q69L"]),
    "parent": avGFP, "scaffold": "avGFP",
    "role": "combined_hero", "expected_tm": 88,
    "notes": "avGFP + sfGFP10 + Q69L (无I152S), 11mut"
})

# --------------------------------------------------------------------
# G1-G3: sfGFP 起点 (htFuncLib 直接)
# 注: sfGFP pos 145=F (Y145F 已有), 所以 Y145M 是 F145M (htFuncLib报告)
# --------------------------------------------------------------------
# G1: sfGFP + I152S + Q69L + S72A (轻量3突变)
candidates.append({
    "name": "G1_sfGFP_I152S_Q69L_S72A",
    "seq": apply(sfGFP, ["I152S", "Q69L", "S72A"]),
    "parent": sfGFP, "scaffold": "sfGFP",
    "role": "combined_balanced", "expected_tm": 90,
    "notes": "sfGFP + I152S + Q69L + S72A 轻量, 3mut"
})

# G2: sfGFP + I152S + Q69L + T108V (轻量3突变 alt)
candidates.append({
    "name": "G2_sfGFP_I152S_Q69L_T108V",
    "seq": apply(sfGFP, ["I152S", "Q69L", "T108V"]),
    "parent": sfGFP, "scaffold": "sfGFP",
    "role": "combined_balanced", "expected_tm": 88,
    "notes": "sfGFP + I152S + Q69L + T108V 轻量 alt, 3mut"
})

# G3: sfGFP + I152S + Q69L + S72A + T108V (中量4突变)
candidates.append({
    "name": "G3_sfGFP_I152S_acid3",
    "seq": apply(sfGFP, ["I152S", "Q69L", "S72A", "T108V"]),
    "parent": sfGFP, "scaffold": "sfGFP",
    "role": "combined_balanced", "expected_tm": 92,
    "notes": "sfGFP + I152S + 完整htFuncLib(Q69L/S72A/T108V), 4mut"
})

# G4: sfGFP + I152S + K166V + Q69L (新组合)
candidates.append({
    "name": "G4_sfGFP_I152S_K166V_Q69L",
    "seq": apply(sfGFP, ["I152S", "K166V", "Q69L"]),
    "parent": sfGFP, "scaffold": "sfGFP",
    "role": "combined_balanced", "expected_tm": 86,
    "notes": "sfGFP + I152S + K166V + Q69L cgre+htFuncLib, 3mut"
})

# --------------------------------------------------------------------
# 验证
# --------------------------------------------------------------------
excl_seqs = set(pd.read_csv(ROOT / "Exclusion_List.csv")["Sequence"].astype(str).str.strip())

print("\n" + "=" * 95)
print(f"{'name':<32} {'scaffold':<8} {'mut':>3} {'chromo':>6} {'len':>4} {'excl':>5} {'status'}")
print("=" * 95)

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
    print(f"{c['name']:<32} {c['scaffold']:<8} {nm:>3} {ch:>6} {len(seq):>4} {excl_mark:>5} {status}")

    c["length"] = len(seq)
    c["n_muts"] = nm
    c["chromophore"] = ch
    c["valid"] = (not issues) and (not in_excl)
    if c["valid"]:
        valid_new.append(c)

print(f"\n通过: {len(valid_new)}/{len(candidates)}")

# 合并候选池
with open(OUT / "candidates_round4_v2.json", encoding="utf-8") as f:
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

with open(OUT / "candidates_round4_v3.json", "w", encoding="utf-8") as f:
    json.dump(existing, f, indent=2, ensure_ascii=False)

print(f"\n保存到 candidates_round4_v3.json (共 {len(existing)} 条)")
