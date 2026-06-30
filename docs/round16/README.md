# Round 16 文档总览

> **日期**: 2026-06-25
> **状态**: ⚠️ 部分完成 — **计算资源探索**
> **核心目标**: 探索更深的 recycles 与多温度采样,以榨干本地算力

---

## 一句话总结

**R16 在 RTX 5080 16GB 上完成的最大规模单次管线扫描**:用 16 重 recycles 重新评估 R15 Top 6 候选,验证项目瓶颈已从"评分函数"转移到"本地算力"。

---

## 文档结构

| 文件 | 用途 |
|------|------|
| [round16_report.md](round16_report.md) | 完整实验报告 |
| [round16_handoff.md](round16_handoff.md) | 交接文档 |
| [round16_runbook.md](round16_runbook.md) | 复现 runbook |
| [round16_pitfalls.md](round16_pitfalls.md) | 踩坑记录 |
| [round16_next_steps.md](round16_next_steps.md) | 下一轮方向 |

---

## 核心数据

| 项目 | 值 |
|:-----|:--|
| 计算设备 | RTX 5080 Laptop 16GB VRAM |
| 候选数 | R15 Top 6 |
| 评估轮数 | 多 recycles (r=8/12/16) |
| 评估时长 | ~6 小时 |
| 最佳排序分 | 0.892(R14 持平) |

---

## 主要发现

1. **r=12 已接近 ESMFold 收敛上限**:继续增加 recycles 没有显著提升
2. **R15 评分体系迁移到 RTX 5080 已失真**:原本的"综合分 0.98"实际上是 fair-esm vs HuggingFace 评分尺度不一致的结果
3. **瓶颈在 GPU 显存**:RTX 5080 16GB 在 r=16 时 batch_size 只能为 1,严重限制吞吐量
4. **关键提示**:要突破 R15 必须切换到云端 A800 80GB

---

*负责人: Trae AI Agent (Claude) | 最后更新: 2026-06-28*
