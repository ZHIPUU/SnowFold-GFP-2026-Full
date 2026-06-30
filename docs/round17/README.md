# Round 17 文档总览

> **日期**: 2026-06-28
> **状态**: ✅ 完成 — **评分体系校准 + 真 Top 1 重定位** 🎯
> **核心成就**: 在 A800 80GB 上用 HuggingFace ESMFold 校准评分,**真实 sort_score 从 0.892 提升到 0.908**
> **关键发现**: **R15 的"综合分 0.98"是评分 bug,真实最佳一直是 0.908**

---

## 文档结构

| 文件 | 用途 |
|------|------|
| [round17_report.md](round17_report.md) | 完整实验报告 |
| [round17_handoff.md](round17_handoff.md) | 交接文档 |
| [round17_runbook.md](round17_runbook.md) | 复现 runbook |
| [round17_pitfalls.md](round17_pitfalls.md) | 踩坑记录 |
| [round17_next_steps.md](round17_next_steps.md) | 下一轮方向 |

---

## 一句话总结

**R17 是项目的"真相回归"轮**:用 A800 80GB 显存 + 正确的 pLDDT 评分尺度(HF 0-1),重新评估 R15 Top 6 候选,得到**真实 Top 1 = R14_A_T01_013 (sort_score=0.908)**。同时验证了 avGFP 路线(改造骨架)的失败。

---

## 核心数据

| 项目 | 值 |
|:-----|:--|
| 计算设备 | **A800 80GB**(首启用) |
| 网络环境 | autodl 网络加速 (ghproxy.com) |
| 模型 | HuggingFace `facebook/esmfold_v1` |
| 评分公式 | `0.4*pTM + 0.3*pLDDT + 0.3*chromo` (pLDDT 0-1) |
| Recycles | 12 |
| 评估时长 | ~5 分钟(R15 Top 6 + 5 新方向) |
| **真实 Top 1** | **R14_A_T01_013** sort_score = **0.908** 🏆 |

---

## 探索方向

| 方向 | 名称 | sort_score | 通过生存底线? |
|:-----|:-----|:---------:|:------------:|
| 1 | sfGFP + S30R 单点 | < 0.60 | ❌ |
| 2 | sfGFP + C端 GIDY | < 0.60 | ❌ |
| 3 | avGFP + sfGFP 关键稳定突变 | < 0.60 | ❌ |
| 4 | sfGFP WT | < 0.60 | ❌ |
| 5 | avGFP WT | < 0.60 | ❌ |

**关键教训**: 任何"非 sfGFP R14 衍生"的序列,在 ESMFold 评分下都无法通过生存底线。

---

## 与历届最佳对比

| 轮次 | 评分尺度 | 真实 sort_score | 虚高显示 |
|:----:|:--------|:--------------:|:-------:|
| R14 | fair-esm 0-100 | 0.892 | 0.892 |
| R15 | fair-esm 0-100 (bug) | 0.892 | **0.983** |
| **R17** | HF 0-1 (校准) | **0.908** | 0.908 |

---

## 主要发现

1. **R14_A_T01_013 是真正的 Top 1**(被 R15 综合分掩盖)
2. **avGFP 路线再次失败**:即使人工引入 sfGFP 6 关键突变,ESMFold 预测 pLDDT 仅 ~0.50
3. **wt 序列作为 MPNN 起点不合适**:sfGFP WT 衍生候选全部 < 0.70
4. **pLDDT 尺度 bug 是项目最大隐患**:R15 之前所有"综合分"对比都有歧义

---

*负责人: Trae AI Agent (Claude) | 最后更新: 2026-06-28*
