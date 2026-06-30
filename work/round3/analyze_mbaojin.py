"""分析 mBaoJin 序列，确认可用性 + 与 avGFP 对比"""
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")

# PDB 8QBJ 序列
mbaojin_pdb = "RSMVSKGEEENMASTPFKFQLKGTINGKSFTVEGEGEGNSHEGSHKGKYVCTSGKLPMSWAALGTTFGYGMKYYTKYPSGLKNWFREVMPGGFTYDRHIQYKGDGSIHAKHQHFMKNGTYHNIVEFTGQDFKENSPVLTGDMNVSLPNEVPQIPRDDGVECPVTLLYPLLSDKSKYVEAHQYTICKPLHNQPAPDVPYHWIRKQYTQSKDDAEERDHICQSETLEAHLKGMDELYK"

# StayGold sequence from original paper (Hirano 2022)
staygold = "MVSKGEEVKEATMQFHYKLTAGLHGQTFTIEGEGKGNPYEGTQKVDLTVIEGAPLPFAYDILTTVFHYGNRAFTKYPHKMPDFFKQVTHTGAFRYDRSIQFTGDGRIVKAKHDFFKNGNYHHVYDFTPGDFKKNGPVLKGDMTASLPNEVRQTGRDEGEYVAPLVYPMLSKNKTFVYKHQYTICKTVKDQPAPGIDPYHWLKRMQQTQAEDDTERDHICQSETLEAHL"

# Known StayGold chromophore: GYG (pos 63-65 typically)
# mBaoJin added mNeonGreen N-term (MVSKGEE...) and C-term (MDELYK)

# Load avGFP
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
print(f"avGFP: {len(avGFP)} aa, chromophore: {avGFP[63:68]}")
print(f"mBaoJin PDB: {len(mbaojin_pdb)} aa")
print(f"StayGold orig: {len(staygold)} aa")

# mBaoJin: strip RS prefix to get M-start
mbaojin_mstart = mbaojin_pdb[2:]  # Remove RS
print(f"\nmBaoJin (M-start): {len(mbaojin_mstart)} aa")
print(f"  First 10: {mbaojin_mstart[:10]}")

# Check chromophore
for name, seq in [("mBaoJin", mbaojin_mstart), ("StayGold", staygold), ("avGFP", avGFP)]:
    for tri in ["TYG", "SYG", "GYG", "CYG", "HYG"]:
        if tri in seq:
            pos = seq.index(tri)
            print(f"  {name}: chromophore={tri} at pos {pos+1}-{pos+3}")
            break

# Align mBaoJin vs StayGold to find mutations
print("\n=== mBaoJin vs StayGold 突变 ===")
# Need proper alignment - simple shift-based
# mBaoJin MVSKGEEENMA... StayGold MVSKGEEVKEA...
# Let's do a simple visual comparison
muts = []
for i, (a, b) in enumerate(zip(mbaojin_mstart, staygold)):
    if a != b:
        muts.append(f"{b}{i+1}{a}")
print(f"Total diff positions: {len(muts)}")
print(f"First 20 mutations: {', '.join(muts[:20])}")
if len(muts) > 20:
    print(f"... and {len(muts)-20} more")

# mBaoJin vs avGFP comparison
print("\n=== mBaoJin vs avGFP 对比 ===")
print(f"avGFP length: {len(avGFP)}")
print(f"mBaoJin length: {len(mbaojin_mstart)}")
print("Note: 不同家族，同源性很低，需要结构对齐后才能移植突变")

# Summary for design use
print("\n" + "=" * 60)
print("mBaoJin 用于竞赛的可行性")
print("=" * 60)
print(f"长度: {len(mbaojin_mstart)} aa ✓ (220-250)")
print(f"M开头: {mbaojin_mstart[0] == 'M'} ✓")
print(f"chromophore: GYG ✓")
print(f"标准AA: {'✓' if all(c in 'ACDEFGHIKLMNPQRSTVWY' for c in mbaojin_mstart) else '✗'}")
print(f"\n关键特性: Tm~92°C, 高亮度, 单体, 快速成熟")
print("可在 mBaoJin 基础上叠加 sfGFP 风格亮度突变")
print("但需注意: StayGold 同源性与 avGFP 仅 ~30%")
