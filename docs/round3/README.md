# Round 3 项目文档 — GFP 蛋白设计

> **项目目标**:为 Synbio Challenges 2026 设计 6 条兼具高 CFPS 亮度和 72°C 热稳定性的 GFP 变体。
>
> **当前状态**:Round 3 已完成,提交 [../submission_yourteamname.csv](../submission_yourteamname.csv) 已就位,6 条候选全部通过排除列表和合规检查,均经过 ESMFold pLDDT 结构验证。
>
> **本文档面向接手者**:请按顺序阅读,理解 Round 3 的改进逻辑与遗留问题。

---

## 快速导航

| 文档 | 内容 | 优先级 |
|------|------|--------|
| [round3_01_overview.md](round3_01_overview.md) | Round 3 总览、相对 Round 1/2 的改进 | 必读 |
| [round3_02_methodology.md](round3_02_methodology.md) | 完整方法论与实现路径 | 必读 |
| [round3_03_results.md](round3_03_results.md) | 成果数据、6 条候选详细解释 | 必读 |
| [round3_04_challenges.md](round3_04_challenges.md) | 难点与坑 — 重点看 ESMFold 部分 | 必读 |
| [round3_05_open_questions.md](round3_05_open_questions.md) | 待解疑点 — 哪些还没解决 | 推荐 |
| [round3_06_next_steps.md](round3_06_next_steps.md) | 下一步方向 — 接手后做什么 | 推荐 |
| [round3_07_handoff.md](round3_07_handoff.md) | 给下一个 AI 的快速上手 | 紧急接手必读 |

---

## TL;DR

**Round 3 关键决策**:
1. ❌ 删除了 Round 2 的 Seq 5 (ppluGFP WT) — **命中排除列表**
2. ✅ 新增 mBaoJin 母体候选 — Tm~92°C, 高亮度, 单体
3. ✅ **首次完成 ESMFold pLDDT 结构验证** — Round 1/2 完全缺失
4. ✅ 全部候选突变数 ≤ 10 — 遵循 Arcadia 2025 论文发现(汉明距离>4 功能急剧下降)
5. ✅ 设计原则从"数据驱动"完全转为"论文知识+结构验证双驱动"

**Round 3 提交 vs Round 2**:
- 6/6 通过全量排除列表(135,414 条) — Round 2 仅 5/6
- 6/6 通过 ESMFold 折叠验证 — Round 2 完全无验证
- 最大突变数 10(Round 3) vs 19(Round 2) — 风险大幅降低

---

## 历史轮次索引

- **Round 1**: 加性 Ridge + ESM2-150M, 数据驱动, 综合分 9.50(乐观) — 见 [../01_achievements.md](../01_achievements.md) §Round 1
- **Round 2**: XGBoost GPU + ESM2-650M + Epistasis, val R²=0.9162 但 OOD 失效 — 见 [../01_achievements.md](../01_achievements.md) §Round 2 和 [round3_03_results.md](round3_03_results.md)
- **Round 3**(当前): 文献调研驱动 + ESMFold 验证 + 严格保守 — 见本目录所有 round3_*.md 文档

---

## 当前 AI session

- **session**: mvs_28b21be1eb074dfa9bcbf4c28733a5fa
- **Root session**: mvs_7706ec6fc7c74872a207498ef6d551ed
- **比赛**: Synbio Challenges 2026 — Protein Design
- **最后更新**: 2026-06-21(Round 3 提交就位)

---

*如果时间紧迫:先读 [round3_07_handoff.md](round3_07_handoff.md) (5 分钟), 再读 [round3_01_overview.md](round3_01_overview.md) (10 分钟)*