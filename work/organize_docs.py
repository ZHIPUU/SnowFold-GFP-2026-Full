"""整理 docs 目录,按轮次归档"""
import shutil
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design\docs")
ROUND3_DIR = ROOT / "round3"

# Round 3 文档移动
round3_files = [
    "round3_README.md",
    "round3_01_overview.md",
    "round3_02_methodology.md",
    "round3_03_results.md",
    "round3_04_challenges.md",
    "round3_05_open_questions.md",
    "round3_06_next_steps.md",
    "round3_07_handoff.md",
]

print("=== 移动 Round 3 文档 ===")
for fname in round3_files:
    src = ROOT / fname
    dst = ROUND3_DIR / fname
    if src.exists():
        shutil.move(str(src), str(dst))
        print(f"  ✓ {fname} → round3/")
    else:
        print(f"  ✗ {fname} 不存在")

# 列出剩余的混合文档(Round 1+2)
print("\n=== 保留在 docs/ 根目录(Round 1+2 混合) ===")
for f in sorted(ROOT.glob("*.md")):
    print(f"  - {f.name}")

# 列出 round3 子目录
print("\n=== docs/round3/ 子目录 ===")
for f in sorted(ROUND3_DIR.glob("*.md")):
    print(f"  - {f.name}")

print("\n完成!")