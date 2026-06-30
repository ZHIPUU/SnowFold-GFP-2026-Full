# Round 12 文档总览

> **日期**: 2026-06-23
> **状态**: ✅ 完成 — **极限 recycles 冲分**
> **核心成就**: 排序分达到 **0.884** (R12 Top 1)
> **策略**: 对 R11 Top 候选做 num_recycles=4/8/12/16/20 扫描

---

## 文档结构

| 文件 | 用途 |
|------|------|
| [round12_report.md](round12_report.md) | 完整实验报告 |
| [round12_handoff.md](round12_handoff.md) | 交接文档 |
| [round12_runbook.md](round12_runbook.md) | 复现 runbook |
| [round12_pitfalls.md](round12_pitfalls.md) | 踩坑记录 |
| [round12_next_steps.md](round12_next_steps.md) | 下一轮方向 |

---

## 一句话总结

**R12 对 R11 Top 6 候选做 num_recycles 扫描，发现 r=8 已是甜蜜点，r=12-20 增益有限；同时验证了 C 端截断 2aa 的隐藏优化维度。**
