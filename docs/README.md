# 2026 Protein Design — GFP 蛋白设计项目文档

> **项目目标**: 为 Synbio Challenges 2026 设计 6 条兼具高 CFPS 亮度和 72°C 热稳定性的 GFP 变体,产出比赛提交 CSV + 设计思路 PDF + GitHub README。
>
> **当前状态**: **Round 19 已完成,真实 sort_score 0.9321 🏆**(R15 虚标 0.98 实际是评分 bug;详见 [round16-19_summary.md](round16-19_summary.md))
>
> **本文档面向接手者**: 请按轮次顺序阅读,理解整个项目演进

---

## 📁 文档结构

```
docs/
├── README.md                  ← 你正在读(总索引)
│
├── 01-08_*.md                 ← Round 1+2 通用文档(混合)
│   ├── 01_achievements.md     ← 两轮成果汇总
│   ├── 02_methodology.md      ← Round 1+2 pipeline 方法论
│   ├── 03_challenges.md       ← Round 1+2 难点与坑
│   ├── 04_open_questions.md   ← Round 2 时代的待解疑点
│   ├── 05_next_steps.md       ← 最新改进方向
│   ├── 06_paper_kb.md         ← 论文突变知识库
│   ├── 07_handoff.md          ← Round 2 时期的接手指南
│   └── 08_appendix.md         ← 文件清单、命令速查、引用
│
├── round3/                    ← Round 3 专属文档
│   ├── README.md
│   ├── round3_01_overview.md
│   ├── round3_02_methodology.md
│   ├── round3_03_results.md
│   ├── round3_04_challenges.md
│   ├── round3_05_open_questions.md
│   ├── round3_06_next_steps.md
│   └── round3_07_handoff.md
│
├── round4/                    ← Round 4 专属文档
│   ├── README.md
│   ├── round4_design_rationale.md
│   ├── round4_next_steps.md
│   ├── round4_pitfalls.md
│   ├── round4_report.md
│   ├── round4_runbook.md
│   └── round4_tech_reference.md
│
├── round5/                    ← Round 5 专属文档
│   └── round5_handoff.md
│
├── round10-15_summary.md      ← R10-R15 6 轮综合总览
├── round10/                   ← Round 10 专属文档 (范式转变 +13.8%)
├── round11/                   ← Round 11 专属文档 (近原点搜索 +23.1%)
├── round12/                   ← Round 12 专属文档 (极限 recycles +23.6%)
├── round13/                   ← Round 13 专属文档 (多样性补救)
├── round14/                   ← Round 14 专属文档 (4路并行 +24.7%)
├── round15/                   ← Round 15 专属文档 (多模型共识 +37.5% 虚标)
│
├── round16/                   ← Round 16 专属文档 (本地算力瓶颈 + 评分 bug 发现)
├── round17/                   ← Round 17 专属文档 (A800 校准 +27.0%)
├── round18/                   ← Round 18 专属文档 (MPNN 大规模 +29.2%)
├── round19/                   ← Round 19 专属文档 (多样性探索 +30.4% 🏆)
└── round16-19_summary.md      ← R16-R19 4 轮综合总览
```

---

## 🚀 快速导航(按角色)

### 🆕 新接手者(15 分钟快速理解)

1. [round16-19_summary.md](round16-19_summary.md) — R16-R19 总览(5 min)
2. [R19 README](round19/README.md) — 当前最佳 (5 min)
3. 看 [D:\workspace\round19\submission_r19.csv](file:///D:/workspace/round19/submission_r19.csv) 验证当前提交 (1 min)

### 📊 看历史背景(30 分钟)

1. [round10-15_summary.md](round10-15_summary.md) — R10-R15 总览 (10 min)
2. [R15 README](round15/README.md) — **注意: R15 综合分 0.98 是评分 bug** (10 min)
3. [round16-19_summary.md](round16-19_summary.md) — R17 校准真相 (10 min)

### 🔬 做研究改进(1-2 小时)

1. [R17 report](round17/round17_report.md) — 评分 bug 修复 (15 min)
2. [R18 report](round18/round18_report.md) — MPNN 大规模 (15 min)
3. [R19 report](round19/round19_report.md) — 多样性探索 (15 min)
4. [R19 next steps](round19/round19_next_steps.md) — 下一轮方向 (15 min)

---

## 📅 项目演进时间线(完整)

### Round 1(已完成,2026-05)
- **方法**: 加性 Ridge + ESM2-150M,纯数据驱动
- **成果**: 综合分 9.50 估计(乐观)

### Round 2(已完成,2026-06-21 前)
- **方法**: XGBoost GPU + ESM2-650M + Epistasis,论文突变堆砌
- **问题**: 模型对 OOD 失效

### Round 3(已完成,2026-06-21)
- **方法**: 文献调研驱动 + ESMFold pLDDT 验证 + 严格保守
- **改进**: 6/6 通过全量排除列表 / 6/6 通过 ESMFold 验证

### Round 4-9(已完成,2026-06-22)
- **方法**: ProteinMPNN 多次尝试
- **最佳**: R4 MPNN_T01_014 (sort_score=0.715)

### Round 10(已完成,2026-06-23) 🚀 范式转变
- **方法**: 迭代修复管线(5 核心固定 + 完整预测 PDB + 多温度)
- **真实分数**: 0.815 (+13.8% vs R4)

### Round 11(已完成,2026-06-23) 🚀🚀 重大突破
- **方法**: R10 最佳作父代,仅重设计 C 端尾(28 突变)
- **真实分数**: 0.881 (+23.1% vs R4)

### Round 12-14(已完成,2026-06-23~24) 微突破
- **真实分数**: R12=0.884, R13=0.884, R14=0.892

### Round 15(已完成,2026-06-25) ⚠️ **评分 bug**
- **方法**: 多 recycles 投票 + ESM2 似然 + 稳定性验证
- **虚标分数**: 0.983 (BUG)
- **真实分数**: 0.892 (= R14,无实质提升)

### Round 16(已完成,2026-06-25) ⚠️ 算力瓶颈暴露
- **发现**: RTX 5080 16GB 是物理上限
- **发现**: fair-esm vs HF pLDDT 尺度 bug

### Round 17(已完成,2026-06-28) 🎯 **评分校准**
- **方法**: A800 80GB + HuggingFace ESMFold
- **真实分数**: **0.908** (+27.0% vs R4)
- **关键**: 修复 R15 评分 bug,真 Top 1 = R14_A_T01_013

### Round 18(已完成,2026-06-28) 🚀🚀🚀 MPNN 大规模
- **方法**: R17 Top 6 父代 × 3 温度 × 50 候选 = 900 候选
- **真实分数**: **0.9240** (+29.2% vs R4)

### Round 19(已完成,2026-06-28) 🏆🏆🏆 多样性探索
- **方法**: 9 父代 × 5 温度 × 150 候选 = 6750 候选
- **真实分数**: **0.9321** (+30.4% vs R4)
- **WT 路线完全失败**:验证近原点是唯一有效路径

---

## 🎯 真实分数演进 (注意 R15 是 bug!)

| 轮次 | 真实 sort_score | 虚标显示 |
|:----:|:---------------:|:-------:|
| R4 | 0.715 | 0.715 |
| R10 | 0.815 | 0.815 |
| R11 | 0.881 | 0.881 |
| R12 | 0.884 | 0.884 |
| R13 | 0.884 | 0.884 |
| R14 | 0.892 | 0.892 |
| R15 | **0.892** | **0.983** ⚠️ BUG |
| R17 | 0.908 | 0.908 |
| R18 | 0.924 | 0.924 |
| **R19** | **0.9321** 🏆 | 0.9321 |

---

## 🔑 关键判断规则

1. **永远用 HuggingFace `out.plddt[..., 1]`** 评分(0-1 尺度)
2. **R15 之前的综合分对比有歧义** — 请用 R17 校准后的真实分数
3. **近原点搜索是唯一有效路径**: WT 路线天花板 0.79
4. **A800 80GB 是必备**: 本地 16GB 无法突破
5. **ProteinMPNN 必须用 ESMFold 预测 PDB 作输入**,非晶体
6. **MPNN fixed_positions key 必须是 basename**,不是全路径
7. **固定 5 核心 [65,66,67,96,222]** 是最优约束
8. **生色团 pLDDT 0.95 = 评分天花板**: 必须换骨架突破

---

## 📂 工作目录

### 本地工作区
```
D:\workspace\
├── round17\    ← R17 重评结果 (sort_score 0.908)
├── round18\    ← R18 MPNN 大规模 (sort_score 0.924)
└── round19\    ← R19 多样性探索 (sort_score 0.9321 🏆 推荐提交)
```

### 项目工作区
```
d:\生信\2026Protein Design\
├── docs\        ← 项目文档
├── work\        ← 历史工作目录
├── 预选序列\    ← 旧提交文件(R15 虚标,需替换)
└── Exclusion_List.csv ← 排除列表(135,414 条)
```

---

## 📞 联络信息

- **比赛**: Synbio Challenges 2026 — Protein Design
- **当前真实最佳**: R19 sort_score 0.9321
- **A800 服务器**: gssh session `9ca7acb1b94c`
- **最后更新**: 2026-06-28 (R19 完成)

---

*文档结构说明: `docs/01-08_*.md` 是 Round 1+2 时代的"轮次混合"文档。`docs/round10-15_summary.md` 总结 R10-R15 突破。`docs/round16-19_summary.md` 总结 R16-R19 真突破(评分校准 + A800 大规模)。R17 是项目分水岭:R17 之前的"分数"实际是 fair-esm bug,R17 之后才是真实分数。*
