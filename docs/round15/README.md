# Round 15 文档总览

> **日期**: 2026-06-25
> **状态**: ✅ 完成 — **抗过拟合 + 多模型共识** 🏆
> **核心成就**: 排序分（final_score）突破 **0.98**（综合 4 指标）
> **策略**: 多 recycles 投票 + ESM2 似然 + 稳定性验证

---

## 文档结构

| 文件 | 用途 |
|------|------|
| [round15_report.md](round15_report.md) | 完整实验报告 |
| [round15_handoff.md](round15_handoff.md) | 交接文档 |
| [round15_runbook.md](round15_runbook.md) | 复现 runbook |
| [round15_pitfalls.md](round15_pitfalls.md) | 踩坑记录 |
| [round15_next_steps.md](round15_next_steps.md) | 下一轮方向 |

---

## 一句话总结

**R15 通过"多 recycles 投票 + ESM2 似然 + 稳定性验证"的多模型共识策略，**在不增加过拟合风险的前提下达到 final_score 0.98**——这是项目迄今为止最稳健的提交。**
