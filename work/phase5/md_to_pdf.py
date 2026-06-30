"""将设计文档 markdown 转 PDF (用 reportlab)"""
from pathlib import Path
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak, Table, TableStyle
from reportlab.lib import colors
from reportlab.lib.units import cm

WORK = Path(r"D:\生信\2026Protein Design\work")
PHASE5 = WORK / "phase5"

with open(PHASE5 / "design_doc.md", "r", encoding="utf-8") as f:
    md_content = f.read()

# 简易 markdown 解析 (支持 # 标题, **粗体**, - 列表, 段落)
def parse_md(text):
    elements = []
    lines = text.split("\n")
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("# "):
            elements.append(("h1", line[2:]))
        elif line.startswith("## "):
            elements.append(("h2", line[3:]))
        elif line.startswith("### "):
            elements.append(("h3", line[4:]))
        elif line.startswith("|") and i + 1 < len(lines) and lines[i+1].startswith("|---"):
            # 表格
            header = [c.strip() for c in line.strip("|").split("|")]
            i += 2
            rows = []
            while i < len(lines) and lines[i].startswith("|"):
                row = [c.strip() for c in lines[i].strip("|").split("|")]
                rows.append(row)
                i += 1
            elements.append(("table", (header, rows)))
            continue
        elif line.startswith("- "):
            elements.append(("li", line[2:]))
        elif line.strip().startswith("```"):
            # 代码块: 收集到下一个 ```
            i += 1
            code = []
            while i < len(lines) and not lines[i].strip().startswith("```"):
                code.append(lines[i])
                i += 1
            elements.append(("code", "\n".join(code)))
        elif line.strip() == "---":
            elements.append(("hr", ""))
        elif line.strip():
            elements.append(("p", line))
        i += 1
    return elements

elements = parse_md(md_content)

# 生成 PDF
doc = SimpleDocTemplate(
    str(PHASE5 / "design_doc.pdf"),
    pagesize=A4,
    rightMargin=2*cm, leftMargin=2*cm,
    topMargin=2*cm, bottomMargin=2*cm
)

styles = getSampleStyleSheet()
title_style = ParagraphStyle('Title', parent=styles['Heading1'], fontSize=18, spaceAfter=12)
h2_style = ParagraphStyle('H2', parent=styles['Heading2'], fontSize=14, spaceBefore=12, spaceAfter=6, textColor=colors.darkblue)
h3_style = ParagraphStyle('H3', parent=styles['Heading3'], fontSize=12, spaceBefore=8, spaceAfter=4, textColor=colors.darkblue)
p_style = ParagraphStyle('P', parent=styles['Normal'], fontSize=10, spaceAfter=4, leading=14)
li_style = ParagraphStyle('Li', parent=styles['Normal'], fontSize=10, leftIndent=20, spaceAfter=2, leading=13)
code_style = ParagraphStyle('Code', parent=styles['Code'], fontSize=8, leftIndent=10, spaceAfter=4, backColor=colors.lightgrey)

flowables = []
for kind, content in elements:
    if kind == "h1":
        flowables.append(Paragraph(content, title_style))
    elif kind == "h2":
        flowables.append(Paragraph(content, h2_style))
    elif kind == "h3":
        flowables.append(Paragraph(content, h3_style))
    elif kind == "p":
        # 简单 markdown 处理: **粗体**
        text = content
        text = text.replace("**", "")  # 简化, 直接去掉
        text = text.replace("`", "")
        flowables.append(Paragraph(text, p_style))
    elif kind == "li":
        text = content.replace("**", "").replace("`", "")
        flowables.append(Paragraph("• " + text, li_style))
    elif kind == "code":
        flowables.append(Paragraph(content.replace("<", "&lt;").replace(">", "&gt;"), code_style))
    elif kind == "hr":
        flowables.append(Spacer(1, 6))
    elif kind == "table":
        header, rows = content
        data = [header] + rows
        t = Table(data, colWidths=[1.5*cm, 1.5*cm, 1.2*cm, 4*cm, 1.2*cm, 1.0*cm, 1.2*cm, 3.5*cm])
        t.setStyle(TableStyle([
            ('BACKGROUND', (0,0), (-1,0), colors.darkblue),
            ('TEXTCOLOR', (0,0), (-1,0), colors.white),
            ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
            ('FONTSIZE', (0,0), (-1,-1), 8),
            ('ALIGN', (0,0), (-1,-1), 'LEFT'),
            ('VALIGN', (0,0), (-1,-1), 'TOP'),
            ('GRID', (0,0), (-1,-1), 0.25, colors.grey),
            ('ROWBACKGROUNDS', (0,1), (-1,-1), [colors.white, colors.lightgrey]),
        ]))
        flowables.append(t)
        flowables.append(Spacer(1, 6))

doc.build(flowables)
print(f"[OK] PDF saved to {PHASE5 / 'design_doc.pdf'}")