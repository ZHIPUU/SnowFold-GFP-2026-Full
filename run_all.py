"""
一键复现脚本: 从原始数据到最终 6 条
================================================
"""
import subprocess
import sys
from pathlib import Path

WORK = Path(r"D:\生信\2026Protein Design\work")

steps = [
    ("Phase 1.2: 数据探索", WORK / "phase1" / "phase1_2_explore.py"),
    ("Phase 1.3: 加性模型", WORK / "phase1" / "phase1_3_v2.py"),
    ("Phase 1.5: 组合搜索", WORK / "phase1" / "phase1_5_search.py"),
    ("Phase 1.6: 筛选 top-144", WORK / "phase1" / "phase1_6_filter.py"),
    ("Phase 2: ESM2 打分", WORK / "phase2" / "phase2_esm_scoring.py"),
    ("Phase 3: 终选 6 条", WORK / "phase3" / "phase3_finalize.py"),
    ("Phase 5: 设计文档", WORK / "phase5" / "gen_design_doc.py"),
    ("Phase 5: PDF", WORK / "phase5" / "md_to_pdf.py"),
]

print("=" * 60)
print("一键复现: 2026 Protein Design GFP 变体设计")
print("=" * 60)

for name, script in steps:
    print(f"\n>>> {name}")
    print(f"    {script}")
    r = subprocess.run([sys.executable, str(script)], capture_output=True, text=True, timeout=1800)
    print(f"    Exit: {r.returncode}")
    if r.returncode != 0:
        print("    STDERR (last 500):")
        print(r.stderr[-500:])
        print(f"\n!!! Failed at {name}")
        sys.exit(1)
    print("    OK")

print("\n" + "=" * 60)
print("ALL DONE")
print("=" * 60)
print(f"Submission: {WORK / 'submission_yourteamname.csv'}")
print(f"Design PDF: {WORK / 'phase5' / 'design_doc.pdf'}")
print(f"Final 6:   {WORK / 'phase3' / 'final_6_candidates.csv'}")