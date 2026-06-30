"""检查 LigandMPNN 输出 chromophore 质量"""
import re, glob

for task in ["sfGFP_lmpnn_T01", "sfGFP_lmpnn_T03", "avGFP_lmpnn_T01"]:
    fa = f"D:/生信/2026Protein Design/work/round5/lmpnn_output/{task}/seqs/"
    files = glob.glob(fa + "*.fa")
    if not files: continue
    with open(files[0]) as f:
        content = f.read()
    blocks = [b for b in content.strip().split(">") if b.strip()]
    wt = blocks[0].split("\n", 1)[1].strip()
    print(f"\n[{task}]")
    print(f"  WT len={len(wt)}")
    for cb in ["TYG", "SYG", "GYG"]:
        if cb in wt: print(f"  chromo {cb} at pos {wt.index(cb)+1}")
    kept = 0
    total = 0
    for b in blocks[1:11]:
        lines = b.split("\n", 1)
        if len(lines) < 2: continue
        seq = lines[1].strip().replace("\n","")
        total += 1
        if any(cb in seq for cb in ["TYG","SYG","GYG"]):
            kept += 1
        else:
            print(f"  设计#{total}: pos 60-75 = {seq[60:75]}")
    print(f"  chromophore 保留: {kept}/{total}")