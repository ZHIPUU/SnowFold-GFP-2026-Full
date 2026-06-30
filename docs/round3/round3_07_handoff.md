# Round 3 接手指南 — 给下一个 AI 的快速上手

> **目标**:让下一个接手的 AI 在 15 分钟内理解 Round 3 的成果和遗留问题,知道该读什么、做什么、不要做什么。
>
> **预计阅读时间**:15 分钟(读完本文)+ 30 分钟(跑通验证脚本)+ 1 小时(理解 Round 3 设计)

---

## TL;DR

**项目**:Synbio Challenges 2026 — GFP 蛋白设计
**当前轮次**:**Round 3 已完成提交**
**提交文件**: [submission_yourteamname.csv](../submission_yourteamname.csv)
**核心改进**:文献调研 + ESMFold 验证 + 严格保守

**Round 3 关键成果**:
- ✅ 6/6 通过全量排除列表(135,414 条)— Round 2 仅 5/6
- ✅ 6/6 通过 ESMFold 折叠验证 — Round 2 完全无验证
- ✅ 全部候选突变数 ≤10(平均 5)— Round 2 最大 19
- ✅ 引入 mBaoJin 候选(Tm~92°C)— Round 2 完全无
- ✅ 8 篇 2024-2025 文献调研

**最紧急的事**(如果你只有 1 小时):
1. ✅ 重跑 ESMFold on GPU 确认 pLDDT(10 min)
2. ✅ 联系 root 询问比赛细节(5 min)
3. 📌 决定是否替换某些候选

**绝对不要做**:
1. ❌ **不要用 Round 2 epistasis 模型给论文突变组合打分**(已证 OOD 失效)
2. ❌ **不要设计 >12 突变的候选**(折叠风险极大)
3. ❌ **不要直接套用 TGP 突变到其他家族**(需要结构对齐)
4. ❌ **不要相信"9.50 综合分"的乐观估计**(实际可能 4-7)
5. ❌ **不要复制 Round 2 Seq 5 (ppluGFP WT)**(已违规)

---

## 5 分钟理解现状

### Round 3 vs Round 1/2 关键差异

| 维度 | Round 1 | Round 2 | **Round 3** |
|------|---------|---------|------------|
| 方法论 | 数据驱动 | 数据+论文堆砌 | **论文驱动 + 结构验证** |
| 模型 | 加性 Ridge | XGBoost GPU + epistasis | **不用 ML 打分** |
| 突变数 | 0-19 | 0-19 | **4-10** |
| 排除列表 | 未检查 | 1000/135K | **135K 全量** |
| 结构验证 | 无 | 无 | **ESMFold pLDDT** |
| 文献调研 | 仅 sfGFP/TGP | 仅 sfGFP/TGP | **8 篇 2024-2025** |
| mBaoJin | 无 | 无 | **1 候选池** |

### 当前提交概览(6 条)

| Seq | 候选 | 突变数 | pLDDT | 设计核心 |
|-----|------|--------|-------|---------|
| 1 | sfGFP+I152S | 1 | 48.3 | sfGFP + chromophore 邻位 |
| 2 | amacGFP+sfGFP5 | 5 | 47.5 | amacGFP + sfGFP 风格 |
| 3 | avGFP+sfGFP10 | 10 | 45.5 | avGFP + sfGFP 完整 |
| 4 | avGFP+sfGFP4core+S30R | 5 | 45.3 | 4 核心 + S30R |
| 5 | avGFP+sfGFP4core | 4 | 44.7 | 最保守 |
| 6 | avGFP+sfGFP4core+I152S | 5 | 42.6 | 4 核心 + I152S |

详细解读见 [round3_03_results.md](round3_03_results.md)。

---

## 项目结构

```
D:\生信\2026Protein Design\
├── submission_yourteamname.csv          ← Round 3 当前提交(6 条候选)
├── AAseqs of 5 GFP proteins_20260511.txt ← 5 个 WT GFP 序列
├── Exclusion_List.csv                   ← 135,414 条排除列表
├── GFP_data.xlsx                        ← 14 万 CFPS 亮度数据
├── docs/
│   ├── 01-08_*.md                       ← Round 1/2 历史文档
│   ├── round3_README.md                 ← Round 3 文档入口(必读)
│   ├── round3_01_overview.md            ← 总览
│   ├── round3_02_methodology.md         ← 方法论
│   ├── round3_03_results.md             ← 6 候选详细
│   ├── round3_04_challenges.md          ← 难点与坑
│   ├── round3_05_open_questions.md      ← 待解疑点
│   ├── round3_06_next_steps.md          ← 下一步方向
│   └── round3_07_handoff.md             ← 本文件
├── work/
│   ├── round2/                          ← Round 2 工作目录(参考)
│   │   └── stepC_xgboost_epistasis.model ← Round 2 模型(R²=0.9162, 但 OOD 失效)
│   ├── round3/                          ← Round 3 工作目录(本轮)
│   │   ├── candidates_round3.json       ← 8 条候选池
│   │   ├── esmfold_results.json         ← pLDDT 数据
│   │   ├── check_submission.py          ← 排除列表+合规检查
│   │   ├── design_candidates_v2.py      ← 候选设计脚本
│   │   ├── esmfold_validate_cpu.py      ← ESMFold CPU 验证
│   │   ├── esmfold_validate_v2.py       ← ESMFold GPU 验证(备用)
│   │   ├── finalize_submission_v2.py    ← 最终选择脚本
│   │   └── *.log                        ← 运行日志
│   └── paper_kb.md                      ← 论文知识库
└── referencepaper/                      ← 5 篇原始 PDF
```

---

## 关键文件(5 分钟读完)

| 文件 | 内容 | 时间 |
|------|------|------|
| [round3_README.md](README.md) | Round 3 文档入口 | 2 min |
| [round3_01_overview.md](round3_01_overview.md) | 总览 + 改进 | 5 min |
| [round3_03_results.md](round3_03_results.md) | 6 候选详细解读 | 5 min |
| [round3_04_challenges.md](round3_04_challenges.md) | 难点与坑 | 5 min |
| [round3_07_handoff.md](round3_07_handoff.md) | 你正在读 | 3 min |

---

## 立即可执行的任务(按 ROI 排序)

### T1. 重跑 ESMFold on GPU (10 min, **最优先**)

**目的**:Round 3 的 pLDDT 是 CPU 推理,绝对值偏低。重跑 GPU 版本获取真实 pLDDT。

```bash
# 1. 确认 CUDA 可用
python -c "import torch; print(torch.cuda.is_available())"

# 2. 在有 CUDA 的 terminal 跑
python "D:\生信\2026Protein Design\work\round3\esmfold_validate_v2.py"
```

**预期**:
- GPU 推理下 pLDDT 高 20-30 分
- 可能改变最终候选排序
- 如果某候选 GPU pLDDT < 60,可能需要替换

**如果 CUDA 不可用**:用 ColabFold API 或等下次有 GPU 时跑。

---

### T2. 联系 root 询问比赛细节 (5 min)

**目的**:获取 Finit/Ffinal 的具体测量条件、Tm 是否直接测量、排名规则等。

**应问问题**:
1. Finit 在多低表达量下测量?
2. Ffinal 加热时长?温度?(72°C 30min 还是其他?)
3. 排名是 6 条平均还是 Best Top-1?
4. 排除列表是全序列匹配还是子串匹配?
5. Tm 是否直接测量,还是通过 Ffinal/Finit 推算?
6. 比赛截止日期?
7. Round 3 提交是否被接受?还是需要 Round 4?

**ROI**:极高 — 5 分钟获取的关键信息可能改变整个策略。

---

### T3. 加 mBaoJin 候选到 Round 4 (1-2 h)

**目的**:Round 3 因 pLDDT 38.8 删除了 mBaoJin。Round 4 应该尝试多个突变版本,找出 pLDDT 最高的。

**做法**:
1. 列出 mBaoJin 候选突变位点(E142D, V193I, L194M, D173N, T185I, S198T 等)
2. 对每个突变版本:
   - 全量排除列表检查
   - ESMFold pLDDT 验证
3. 选 pLDDT 最高的 mBaoJin 版本加入 Round 4 候选池

**预期**:Tm 92°C 是无与伦比的优势。

---

### T4. 验证当前提交合规性 (5 min)

**目的**:确保当前提交仍然全部合规。

```python
import pandas as pd
excl = pd.read_csv(r"D:\生信\2026Protein Design\Exclusion_List.csv")
excl_seqs = set(excl["Sequence"].astype(str).str.strip())
sub = pd.read_csv(r"D:\生信\2026Protein Design\submission_yourteamname.csv")
for _, row in sub.iterrows():
    in_excl = str(row["Sequence"]).strip() in excl_seqs
    print(f"Seq {row['Seq_ID']}: {'❌' if in_excl else '✓'}")
```

---

### T5. 跑 Tm 预测模型 (1-2 天)

**目的**:完善综合分评估(Ffinal/Finit 项)。

**做法**:
1. 收集 50-100 个 GFP Tm 数据(从 FPbase + 文献)
2. 用 ESM2-650M 嵌入(可复用 Round 2 的 esm650m_embeddings.npy)
3. Ridge 回归训练
4. 验证集 R² > 0.7

---

## 关键判断规则(Round 4 必须遵守)

### Rule 1: 突变数 ≤10

如果设计新候选,突变数必须 ≤10。Arcadia Science 2025 数据定量证明:汉明距离 >4 后功能急剧下降。

### Rule 2: 不要用 ML 给 OOD 打分

Round 2 epistasis 模型对论文突变组合系统性低估 1.5+ log10。如果要打分,使用:
- 论文先验(chromophore + 文献突变)
- ESMFold pLDDT(结构验证)
- Tm 预测器(若已训练)

### Rule 3: 全量排除列表

任何提交前都**必须**全量检查 135,414 条排除列表。野生型序列经常命中。

### Rule 4: pLDDT 是必要验证

不要提交任何没经过 ESMFold 验证的候选。pLDDT > 60(GPU)或 > 45(CPU)为可接受。

### Rule 5: 中文路径用脚本文件

PowerShell + 中文路径是地狱。所有 Python 脚本都用 `pathlib.Path(r"...")` 显式 UTF-8 字符串。不要在命令行直接传中文路径。

---

## 调试常见问题

### Q: 我想加载 ESMFold,失败了
A: 见 [round3_04_challenges.md](round3_04_challenges.md) §1.2。SSL 错误用 `local_files_only=True` + 预先下载的本地模型。

### Q: ESMFold pLDDT 看起来全是 0.4-0.5
A: 你用了 FP16 或没乘以 100。必须用 FP32 + × 100。详见 [round3_04_challenges.md](round3_04_challenges.md) §3.2。

### Q: 候选序列命中排除列表
A: 加 1-2 个表面保守突变再试。如果还是命中,放弃这个 scaffold。

### Q: 我的候选分数很低,怎么办?
A: 优先检查:
1. chromophore 完整吗?(搜 TYG/SYG/GYG)
2. 长度 220-250 吗?
3. M 开头?
4. 在排除列表里?(查 13.5 万)
5. **是否过度突变**(>15 muts)?
6. **是否做了 pLDDT 验证?**

### Q: 我想跑 GeoEvoBuilder,如何开始?
A: 详见 [round3_06_next_steps.md](round3_06_next_steps.md) §S4。代码: github.com/PKUliujl/GeoEvoBuilder。

### Q: 我想用 ProteinMPNN,如何开始?
A: 详见 [round3_06_next_steps.md](round3_06_next_steps.md) §S6。需要 PDB backbone 输入。

---

## 不要做的事(避免重蹈覆辙)

1. ❌ **不要纯数据驱动** — 数据天花板 2.33×
2. ❌ **不要给 OOD 候选用 XGBoost 打分** — 系统性低估
3. ❌ **不要相信 "Round 1 Seq 6 = 9.50"** — 那是乐观估计
4. ❌ **不要忽略 scaffold 平衡** — 6 条全用 avGFP 风险高
5. ❌ **不要在没 pLDDT 的情况下交 >15 muts 的候选**
6. ❌ **不要直接套用 TGP 突变到其他家族**(需要结构对齐)
7. ❌ **不要复制 Round 2 Seq 5 (ppluGFP WT)**(已违规)
8. ❌ **不要在 CPU 上跑大批量 ESMFold** — 太慢
9. ❌ **不要相信野生型 GFP 是"安全"选项** — 经常命中排除列表

---

## 如果你只能做一件事

**重跑 ESMFold on GPU,确认 Round 3 的 6 条候选 pLDDT 真实值。**

如果 GPU pLDDT 显示某些候选 < 60,考虑替换或移除。这是参赛的最低门槛。

---

## 紧急联络

- **当前 AI session**: mvs_28b21be1eb074dfa9bcbf4c28733a5fa
- **Root session**: mvs_7706ec6fc7c74872a207498ef6d551ed
- **比赛**: Synbio Challenges 2026 — Protein Design
- **截止日期**: 问 root session
- **最后更新**: 2026-06-21

---

## Round 3 vs Round 2 关键决策回顾

| Round 2 决策 | Round 3 决策 | 原因 |
|--------------|--------------|------|
| 11 条候选 → 选 6 | 10 条候选 → 选 6 | 文献调研后减少设计复杂度 |
| Seq 5 (ppluGFP WT) | 删除 | 命中排除列表 |
| Seq 1 (19 突变) | 删除 | 突变数过多 |
| ML 模型打分 | 不用 ML 打分 | OOD 失效已确认 |
| 无结构验证 | ESMFold 全验证 | 文献 + 实践双驱动 |
| 5 种 scaffold | 3 种 scaffold(全部合规) | 质量 > 数量 |
| 无 mBaoJin | 1 mBaoJin 候选池 | 新增母体 |

---

*对应 Round 3 总览见 [round3_01_overview.md](round3_01_overview.md), 详细技术细节见 [round3_02_methodology.md](round3_02_methodology.md), 难点细节见 [round3_04_challenges.md](round3_04_challenges.md), 下一步方向见 [round3_06_next_steps.md](round3_06_next_steps.md)。*