# 07 · 给下一个 AI 的快速上手指南(Handoff)

> **目标**:让下一个接手的 AI 在 15 分钟内理解项目状态,知道该读什么、做什么、不要做什么。
>
> **预计阅读时间**:15 分钟(读完本文)+ 30 分钟(跑通验证脚本)+ 1 小时(理解 Round 2 设计)

---

## TL;DR

**项目**:Synbio Challenges 2026 — GFP 蛋白设计。**Round 1 完成**(综合分 9.50 估计),**Round 2 已提交**(4 scaffold 覆盖 + 论文突变组合,无 pLDDT 验证)。

**当前最好的资产**:
- `submission_yourteamname.csv` — Round 2 提交
- `work/round2/stepC_xgboost_epistasis.model` — val R²=0.9162 的 epistasis 模型
- `work/paper_kb.md` — 5 篇论文的突变清单
- 8 个 markdown 文档(就是你正在读的 docs/)

**最紧急的事**(如果你只有 1 小时):
1. ✅ 验证 `submission_yourteamname.csv` 全合规(长度、M 开头、AA)
2. ✅ 全量检查排除列表(目前只查了 1000/135000)
3. 📌 如有 GPU 时间,下载 ESMFold 验证 6 条候选结构

**绝对不要做**:
1. ❌ 不要让 XGBoost/epistasis 模型给"论文突变组合"候选打分(模型对 OOD 失效,会全部低估)
2. ❌ 不要尝试再嵌入 ESM2-3B(除非有 4+ 小时空闲)
3. ❌ 不要重新走"纯数据驱动"路线(天花板 2.33×)

---

## 5 分钟理解现状

### 项目结构
```
D:\生信\2026Protein Design\
├── README.md                       ← 项目入口(必读)
├── submission_yourteamname.csv     ← Round 2 提交
├── AAseqs of 5 GFP proteins_*.txt  ← 5 个 WT 序列
├── Exclusion_List.csv              ← 13.5 万排除列表
├── GFP_data.xlsx                   ← 14 万 CFPS 亮度数据
├── docs/                           ← 8 个文档(必读)
├── work/                           ← 所有工作输出
│   ├── round2/                     ← Round 2 主目录
│   │   ├── stepC_xgboost_epistasis.model  ← 主力模型
│   │   ├── candidates_round2_design.csv   ← 11 条候选
│   │   ├── stepC_summary.json             ← 模型指标
│   │   └── *.py                          ← 所有脚本
│   └── paper_kb.md                 ← 论文突变知识库
└── referencepaper/                 ← 5 篇 PDF
```

### Round 1 vs Round 2 差异
- **Round 1**:加性 Ridge + ESM2-150M,数据驱动,综合分 9.50(乐观)
- **Round 2**:XGBoost GPU + ESM2-650M + epistasis 特征,val R²=0.92 但 OOD 失效,最终回退到论文驱动

### 关键文件(5 分钟读完)
| 文件 | 内容 | 时间 |
|---|---|---|
| `README.md` | 项目入口 | 2 min |
| `docs/01_achievements.md` | 两轮成果 | 5 min |
| `docs/06_paper_kb.md` | 论文突变 | 5 min |
| `docs/07_handoff.md` | 你正在读 | 3 min |

---

## 立即可执行的任务(按 ROI 排序)

### T1. 全量排除列表检查(10 min, 必做)

**目的**:确认 6 条候选不在 13.5 万排除列表中。

```python
import pandas as pd
excl = pd.read_csv(r"D:\生信\2026Protein Design\Exclusion_List.csv")
excl_seqs = set(excl.iloc[:, 0].astype(str))
sub = pd.read_csv(r"D:\生信\2026Protein Design\submission_yourteamname.csv")
for _, row in sub.iterrows():
    in_excl = row["Sequence"] in excl_seqs
    print(f"Seq {row['Seq_ID']}: {'❌ IN EXCLUSION!' if in_excl else '✓ OK'} ({len(row['Sequence'])} aa)")
```

**如果命中**:替换该候选(从 work/round2/candidates_round2_design.csv 里选另一条)。

### T2. 加 2 条保守备选(30 min)

**目的**:降低 Seq 1 (19 muts) 风险。

**做法**:编辑 `work/round2/step3c_design_candidates.py`,加 2 条候选:
- avGFP + sfGFP 完整 11 突变(**仅 11 muts**,无 TGP 额外) — 安全牌
- avGFP + sfGFP 5 折叠核心 + TGP 4 稳定核心(共 **9 muts**) — 轻量版

然后从 13 条候选里重新选 6 条(优先保守 + 高 chromophore score)。

### T3. 跑一遍 Round 2 模型打分(5 min)

**目的**:验证模型对 in-distribution variants 的表现(确认不是 bug)。

```bash
python "D:\生信\2026Protein Design\work\round2\step3a_score_141k.py"
```

**预期输出**:max Finit_rel ≈ 2.33× (avGFP, E5K:D75N:N104Y)。

### T4. 验证嵌入文件可加载(2 min)

```python
import numpy as np
emb = np.load(r"D:\生信\2026Protein Design\work\round2\esm650m_embeddings.npy", mmap_mode="r")
print(emb.shape, emb.dtype)  # 应是 (141144, 1280) float32
```

---

## 关键判断规则

### Rule 1: 模型对 OOD 候选失效

如果你设计了一条"全新突变组合"候选,**不要让 XGBoost/epistasis 模型打分**。它们的预测会系统性低估。

**正确做法**:基于 `docs/06_paper_kb.md` 的先验 + chromophore 完整性 + scaffold 平衡。

### Rule 2: 数据天花板 ~2.33×

不要尝试纯数据驱动路线。**任何从老 top_candidates.csv 出发的工作都是死路**(已验证)。

### Rule 3: Tm 没有预测模型

综合分公式中 Ffinal/Finit 没有预测能力。**估算时用文献已知 Tm**:
- avGFP WT ~64°C
- sfGFP ~78°C
- TGP, StayGold, mBaoJin ~85-92°C
- cgreGFP, ppluGFP ~75-80°C(估计)

### Rule 4: pLDDT 验证是金标准

如果有 ESMFold 权重,**必须跑 6 条候选的 pLDDT**。pLDDT < 70 的候选大概率不折叠。

### Rule 5: 中文路径 + PowerShell 是地狱

所有 Python 脚本都用 `pathlib.Path(r"...")` 显式 UTF-8 字符串。`bash` 工具里调 PowerShell 时,经常编码损坏。

---

## 调试常见问题

### Q: 我想加载 ESMFold,失败了
A: 见 `docs/03_challenges.md` §3。SSL 问题绕不过去。**离线下载 14 GB 是唯一方案**,否则放弃。

### Q: XGBoost 训练太慢
A: 已用 GPU (`device='cuda'`)。如果慢,检查 `tree_method='hist'` 和 batch 模式。

### Q: 嵌入重跑要 80+ 分钟
A: 没办法,141K × 1280-d 必须这么慢。可以减小 batch 到 64,但会慢一些。

### Q: 怎么知道哪些突变是"已知好"的?
A: 查 `docs/06_paper_kb.md`。最关键的 6 个: **S65T, S30R, F64L, F99S, M153T, V163A**。

### Q: 我的候选分数很低,怎么办?
A: 优先检查:
1. chromophore 完整吗?(搜 TYG/SYG/GYG)
2. 长度 220-250 吗?
3. M 开头?
4. 在排除列表里?(查 13.5 万)
5. **是否过度突变**(>15 muts)?

---

## 不要做的事(避免重蹈覆辙)

1. ❌ **不要纯数据驱动** — 数据天花板 2.33×
2. ❌ **不要给 OOD 候选用 XGBoost 打分** — 系统性低估
3. ❌ **不要相信 "Round 1 Seq 6 = 9.50"** — 那是乐观估计,实际可能 4-7
4. ❌ **不要忽略 scaffold 平衡** — 6 条全用 avGFP 风险高
5. ❌ **不要在没 pLDDT 的情况下交 >15 muts 的候选** — 可能折叠失败
6. ❌ **不要尝试 ESMFold 加载** — SSL 解决不了,改用其他结构验证
7. ❌ **不要重复造轮子** — 11 条候选已经在 `candidates_round2_design.csv`,直接用

---

## 如果你只能做一件事

**确保 `submission_yourteamname.csv` 6 条候选都通过基本验证 + 排除列表 + chromophore 完整**。这是参赛的最低门槛,其他都是优化。

---

## 紧急联络

- **当前 AI session**:mvs_28b21be1eb074dfa9bcbf4c28733a5fa
- **Root session**:mvs_7706ec6fc7c74872a207498ef6d551ed
- **比赛**:Synbio Challenges 2026 — Protein Design
- **截止日期**:问 root session

---

*对应项目完整说明见 `README.md`,对应技术细节见 `docs/02_methodology.md`,对应未解疑点见 `docs/04_open_questions.md`*