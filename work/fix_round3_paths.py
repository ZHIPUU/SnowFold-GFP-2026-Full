"""批量修复 round3 子目录内文档的相对路径"""
from pathlib import Path

ROUND3 = Path(r"D:\生信\2026Protein Design\docs\round3")

# 修复项: (old, new) 相对路径
fixes = [
    # 指向 docs/ 根的文档,需要加 ../
    ("(01_achievements.md)", "(../01_achievements.md)"),
    ("(02_methodology.md)", "(../02_methodology.md)"),
    ("(03_challenges.md)", "(../03_challenges.md)"),
    ("(04_open_questions.md)", "(../04_open_questions.md)"),
    ("(05_next_steps.md)", "(../05_next_steps.md)"),
    ("(06_paper_kb.md)", "(../06_paper_kb.md)"),
    ("(07_handoff.md)", "(../07_handoff.md)"),
    ("(08_appendix.md)", "(../08_appendix.md)"),
    # round3_README.md → README.md(同名在 round3/ 内)
    ("(round3_README.md)", "(README.md)"),
    # work/ 路径在 docs/round3/ 内需要 ../../work/
    ("(../../work/round3/", "(../../work/round3/"),  # 已经是正确的
]

# 实际 work/ 路径应该是 ../../work/
extra_fixes = [
    ("(work/round3/", "(../../work/round3/"),
]

# 修正 docs/round3/ 子目录中的 work/ 引用
for f in sorted(ROUND3.glob("*.md")):
    text = f.read_text(encoding="utf-8")
    original = text

    # 应用基础修复
    for old, new in fixes:
        text = text.replace(old, new)

    # 应用 work/ 路径修复
    for old, new in extra_fixes:
        text = text.replace(old, new)

    if text != original:
        f.write_text(text, encoding="utf-8")
        print(f"  ✓ 修复 {f.name}")
    else:
        print(f"  - {f.name} 无需修改")

print("\n完成!")