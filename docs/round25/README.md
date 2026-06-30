# Round 25 文档总览

> **日期**: 2026-06-30
> **状态**: ✅ 完成 — **sort_score = 0.9477 🏆 项目历史最高**
> **突破**: pTM/chromo/pLDDT 三项全部创新高

---

## 📁 文档结构

| 文件 | 用途 |
|------|------|
| [round25_report.md](round25_report.md) | 完整实验报告 |
| [round25_handoff.md](round25_handoff.md) | 交接文档 |
| [round25_runbook.md](round25_runbook.md) | 复现 runbook |
| [round25_pitfalls.md](round25_pitfalls.md) | 踩坑记录 |
| [round25_next_steps.md](round25_next_steps.md) | 下一轮方向 |

---

## 📊 核心数据

| 项目 | 值 |
|:-----|:--|
| 父代 | R24 Top 6 (sort_score 0.9447) |
| 温度 | [0.05, 0.1, 0.2, 0.4, 0.6] |
| 总候选 | 6000 |
| Recycles | 8 (screen) + 20 (precise) |
| **Top 1 sort_score** | **0.9477** 🏆 |
| **Top 1 pTM** | **0.9321** |
| **Top 1 chromo** | **0.970** |
| 总耗时 | 12.4 小时 |
| 自动启动 | watcher 脚本自动接力 |

### R25 Top 6

| Seq | Score | pTM | pLDDT | Chromo | Parent |
|:---:|:-----:|:---:|:-----:|:------:|:-------|
| 1 | **0.9477** | **0.9321** | **0.947** | **0.970** | r25_p4 |
| 2 | 0.9468 | 0.9320 | 0.946 | 0.968 | r25_p5 |
| 3 | 0.9468 | 0.9312 | 0.947 | 0.968 | r25_p4 |
| 4 | 0.9466 | 0.9309 | 0.946 | 0.968 | r25_p5 |
| 5 | 0.9464 | 0.9290 | 0.947 | 0.969 | r25_p4 |
| 6 | 0.9463 | 0.9301 | 0.949 | 0.966 | r25_p6 |

---

## 📂 文件

- `D:\workspace\round25\submission_r25.csv` — **推荐提交**
- `D:\workspace\round25\final_6_r25.json` — Top 6 详细
- `D:\workspace\round25\all_passed.json` — 全部通过候选

---

*负责人: Trae AI Agent (Claude) | 完成时间: 2026-06-30 20:42*
