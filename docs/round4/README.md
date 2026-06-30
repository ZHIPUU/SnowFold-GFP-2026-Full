# Round 4 文档目录 (接手者必读)

> **生成时间**: 2026-06-22
> **状态**: 完整 (6/6 文档)
> **适用**: Round 5 接手者 / Round 4 复现 / 设计 PDF 写作

---

## 📑 文档列表 (按阅读顺序)

### 1. [round4_report.md](round4_report.md) - **主报告 (30 分钟读完)**
- 完整 TL;DR: 最终 Top-6 + 三场景预测
- 相对 Round 3 的五大改进
- 三轮 R² 演变
- 完整文件清单 + 关键产出

### 2. [round4_design_rationale.md](round4_design_rationale.md) - **设计思路 (30 分钟读完)**
- 6 候选的"为什么是它"
- 关键设计选择的"为什么"
- 与比赛规则的映射
- 给 Round 5 的设计建议

### 3. [round4_tech_reference.md](round4_tech_reference.md) - **技术参考 (按需查)**
- 硬件 / 软件 / 模型规格
- ESMFold / ProteinMPNN 详细使用
- 5 个 GFP PDB 信息
- 关键文献速查

### 4. [round4_pitfalls.md](round4_pitfalls.md) - **踩坑指南 (20 分钟读完)**
- 15 个踩坑 + 解决方案
- 中文路径 / CPU torch / X 占位符等
- 调试模板 + 关键文件清单

### 5. [round4_runbook.md](round4_runbook.md) - **运行指南 (复现用)**
- 从零复现 Round 4 的 step-by-step
- 完整时间表 (3.3 小时)
- 验证清单 + 紧急恢复

### 6. [round4_next_steps.md](round4_next_steps.md) - **下一步战略 (20 分钟读完)**
- ROI 排序 + 决策树
- 立即必做 / 短期 / 中期 / 长期
- 风险评估 + 给 Round 5 的核心建议

---

## 🚀 快速开始 (15 分钟摘要)

### 最终提交 (推荐)
📁 `D:\生信\2026Protein Design\work\round4\submission_round4_v5.csv`

### Top-6 候选
| Seq | 候选 | 骨架 | pLDDT | 预测综合分 |
|-----|------|------|-------|-----------|
| 1 | MPNN_T01_014 | sfGFP_MPNN | 68.3 ⭐ | 1.04 |
| 2 | MPNN_av_T03_v2_001 | avGFP_MPNN | 61.5 | 0.81 |
| 3 | G1_sfGFP_I152S_Q69L_S72A | sfGFP | 48.6 | 0.90 |
| 4 | H1_avGFP_sfGFP_acid3_I152S | avGFP | 49.7 | 🏆 1.23 |
| 5 | Z3_amacGFP_sfGFP5_I152S | amacGFP | 45.5 | 0.50 |
| 6 | M7_mBaoJin_K173R_Y196F | mBaoJin | 39.0 | 0.54 |

### 三场景预测 (Best Top-1)
- 🔴 悲观: 0.68
- 🟡 中性: **1.23**
- 🟢 乐观: 1.90

### 必做 (未完成项)
1. **设计思路 PDF** (1-2 小时) - 比赛必交
2. **GitHub 仓库** (30 分钟) - 比赛必交
3. **联系 root 询问比赛细节** (5 分钟) - 高 ROI

### 不要做
1. ❌ 用 ML 打分 OOD (Round 2 教训)
2. ❌ 在中文路径跑 ProteinMPNN
3. ❌ 改 v5 提交 (已最优)
4. ❌ 装 CPU 版 PyTorch
5. ❌ 追求 100% 完美 (BP Top-1 已足够)

---

## 🔑 关键发现 (10 条)

1. **ProteinMPNN 在 GFP 上 pLDDT 可达 68** - 显著超过所有手工设计
2. **MPNN_T01_014 是全场最稳候选** (pLDDT 68, pTM 0.765)
3. **mBaoJin Tm 92°C 是最佳热稳奖** - 但 pLDDT 39 风险
4. **5 骨架多样性显著降低风险** - 优于 Round 3 的 5/6 几乎相同
5. **Tm 估值按 pLDDT 分级** - v5 改进让预测更准
6. **H1 (avGFP+13mut) 综合分 1.23 最高** - 但 pLDDT 49 略低
7. **Round 2 ML 失效是领域级陷阱** - 不要用 ML 打分 OOD
8. **CPU vs GPU pLDDT 几乎一致** - 装 torch 时要用对版本
9. **中文路径让 ProteinMPNN 失败** - 必须在 C:\Temp\ 跑
10. **ProteinMPNN 输出含 X 占位符** - 需要后处理填回

---

## 📞 联系信息

### 比赛必交材料 (Round 4 状态)
- ✅ 序列文件: `work/round4/submission_round4_v5.csv`
- ❌ 设计思路 PDF: **待写** (1-2 小时)
- ❌ GitHub 仓库: **待建** (30 分钟)

### AI Session 信息
- Root session: mvs_7706ec6fc7c74872a207498ef6d551ed
- Current session: mvs_28b21be1eb074dfa9bcbf4c28733a5fa
- 比赛: SynBio Challenges 2026 — Protein Design

### 关键参考
- 比赛规则 PDF: `2026Protein Design in Synbio challenges.pdf`
- 5 个 WT 序列: `AAseqs of 5 GFP proteins_20260511.txt`
- 官方数据: `GFP_data.xlsx` (14万条 CFPS)
- 排除列表: `Exclusion_List.csv` (13.5万条)

---

## 📊 文档统计

| 文档 | 字数 | 阅读时间 | 类型 |
|------|------|----------|------|
| round4_report.md | ~6000 | 30 min | 主报告 |
| round4_design_rationale.md | ~5000 | 30 min | 设计思路 |
| round4_tech_reference.md | ~4000 | 按需 | 技术参考 |
| round4_pitfalls.md | ~4000 | 20 min | 踩坑指南 |
| round4_runbook.md | ~3000 | 复现时 | 运行指南 |
| round4_next_steps.md | ~4000 | 20 min | 战略指南 |
| **合计** | **~26000 字** | **~2 小时** | |

---

**完成日期**: 2026-06-22
**作者**: Trae AI Agent (Claude Sonnet)
**目的**: 让下一个 AI 在最短时间内掌握 Round 4 全貌
