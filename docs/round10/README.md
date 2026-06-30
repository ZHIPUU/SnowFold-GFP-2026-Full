# Round 10 文档总览

> **日期**: 2026-06-23
> **状态**: ✅ 完成 — 实现历史性突破（排序分 +13.8% vs R4）
> **核心成就**: 排序分从 0.715 → **0.815**（+0.100, +13.8%）
> **策略**: 迭代式 pLDDT 靶向修复管线（ProteinMPNN + ESMFold 闭环）

---

## 文档结构

| 文件 | 用途 |
|------|------|
| [round10_report.md](round10_report.md) | 完整实验报告（方法、数据、结论） |
| [round10_runbook.md](round10_runbook.md) | 复现 runbook（一步步操作指南） |
| [round10_pitfalls.md](round10_pitfalls.md) | 踩坑记录（给后续 Agent 的警告） |
| [round10_design_rationale.md](round10_design_rationale.md) | 设计哲学与策略选择理由 |
| [round10_next_steps.md](round10_next_steps.md) | 下一轮探索方向 |
| [round10_tech_reference.md](round10_tech_reference.md) | 技术参考（参数、性能、API） |
| [round10_handoff.md](round10_handoff.md) | 交接文档（核心总结） |

---

## 一句话总结

**R10 通过"5 核心残基固定 + 完整预测 PDB 骨架 + 多温度采样"的迭代修复管线，超越了 R4 保持 6 轮的最佳成绩（MPNN_T01_014, pLDDT=68.3）。**
