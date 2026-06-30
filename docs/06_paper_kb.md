# 06 · 论文突变知识库 (Paper Knowledge Base)

> **最新更新**: 2026-06-28 (R19 完成 + 2025 新文献整合)
> 汇总 GFP 相关经典论文 + 2025 新进展,用于补强评分体系和设计策略

---

## 一、GFP 经典文献(保留)

### 1. Superfolder GFP (sfGFP) — Pédelacq et al. 2006
- avGFP + 11 突变 → Tm ~78°C (vs avGFP 64°C)
- **关键突变**: S30R (+1.25 kcal/mol), F64L, S65T, F99S, M153T, V163A, Y39N, N105T, Y145F, I171V, A206V
- **启示**: S65T 是几乎所有现代 GFP 的标准突变

### 2. TGP — Close et al. 2015 (Los Alamos)
- 7 个稳定突变 + F99S/Y145F 核心
- Tm 92°C (热稳定冠军)

### 3. mCherry/mEGFP 系列
- mCherry 起源来自 DsRed 突变
- 亮度/成熟度优化

### 4. sfGFP 在 CFPS 系统
- sfGFP 在 Cell-Free Protein Synthesis 系统中表达良好
- **参赛项目验证场景**: CFPS 表达 + 72°C 热处理

---

## 二、2025 最新文献(本次新增)

### A. RFdiffusion3 (RFD3) — Baker Lab 2025.09

**引用**: [De novo Design of All-atom Biomolecular Interactions with RFdiffusion3](https://pmc.ncbi.nlm.nih.gov/articles/PMC12458353/)

**核心**:
- 全原子扩散模型(包含配体/核酸/小分子)
- 在 RoseTTAFold 基础上扩展
- **首次实现 GFP β-barrel 的 de novo 设计**
- 包含 de novo 跨膜纳米孔 (200-500 pS 导电性)

**对本项目的应用**:
```bash
# 完整生成全新 β-barrel 骨架
git clone https://github.com/RosettaCommons/foundry
# 用 RFD3 生成 10-20 个新 β-barrel
# 然后用 ProteinMPNN 重新设计序列
# 用 ESMFold 验证
```

**预期收益**: 突破 sfGFP 骨架的限制,生色团 pLDDT 可能突破 0.95

### B. Parametric β-barrel Design — Baker Lab PNAS 2025.09

**引用**: [DOI 10.1073/pnas.2425459122](https://pnas.org/doi/10.1073/pnas.2425459122)

**核心**:
- 用参数化圆柱体 + RFdiffusion 引导生成 β-barrel
- 在多个 β-sheet 参数化范围内成功
- X-ray 晶体结构验证 (atomic accuracy)
- 跨膜纳米孔设计 (200-500 pS)

**对本项目的应用**:
- 用参数化方法生成"sfGFP-like"的圆柱骨架
- 在骨架基础上设计生色团

### C. ESM3 GFP "5 亿年进化" — EvolutionaryScale Science 2025

**引用**: [Simulating 500 million years of evolution with a language model](https://www.science.org/doi/10.1126/science.ads0018)

**核心**:
- ESM3 (7B 参数) 多模态生成语言模型
- 直接以生色团 6 关键残基 (Thr62/Thr65/Tyr66/Gly67/Arg96/Glu222) 为提示
- 生成 GFP 与已知 GFP 序列同一性仅 **58%**
- **"等价于 5 亿年进化"**
- 合成 GFP 在体外实验验证保留荧光功能

**对本项目的应用**:
```python
# 直接以生色团为约束生成全新 GFP 序列
from esm.models.esm3 import ESM3
model = ESM3.from_pretrained("esm3-medium-2024-08")

# 关键提示: 6 个生色团残基固定
prompt = ["", "", "", "T", "T", "Y", "G", "", ..., "R", "", ..., "E", "", ""]
generated = model.generate(prompt)
```

**预期收益**: 可能是项目最大的突破机会(完全跳出 sfGFP 限制)

### D. ThermoMPNN ΔΔG 预测 — Kuhlman Lab PNAS 2024

**引用**: [Transfer learning to leverage larger datasets for improved prediction of protein stability changes](https://pmc.ncbi.nlm.nih.gov/articles/PMC10402116/)

**核心**:
- 在 ProteinMPNN 基础上迁移学习
- Megascale 数据集: 180 万测量,300+ 蛋白的每点突变
- 输出 ΔΔG (kcal/mol) — 越负越稳定
- 推理速度: 秒级完成全 SSM (单点突变扫描)

**对本项目的应用**:
```bash
pip install git+https://github.com/Kuhlman-Lab/ThermoMPNN.git
```

```python
from thermompnn import ThermoMPNN
model = ThermoMPNN.from_pretrained()
# 对 R19 Top 6 评估 ΔΔG
ddg = model.predict(pdb_path="/root/autodl-tmp/r19/pdbs/R19_Top1.pdb")
# 选择 ΔΔG < 0 (稳定化) 的候选
```

**预期收益**: 直接填补 Tm 预测空白,综合分 Tm 维度提升

### E. Arcadia GFP CNN Ensemble — 2025.09

**引用**: [Efficient GFP variant design with a simple neural network ensemble](https://thestacks.org/publications/result-gfp-variant-design-nn)

**核心**:
- ESM-2 (15B 参数) embedding + 简单 CNN ensemble
- 在 Arcadia DMS 数据集训练
- 生成 1000 条新 GFP,选 10 条验证
- 实验确认多个候选亮于 baseline
- **完整开源**: [GitHub 仓库](https://github.com/Arcadia-Science/2025-GFP-variant-design) + [Zenodo 数据](https://zenodo.org/records/17088257)
- **完整 RFP-GFP 双报告系统**实验流程

**对本项目的应用**:
- 直接用 Arcadia 训练好的 CNN 模型给候选打分
- 4 个月完成从概念到实验验证的全流程参考

### F. Steering Protein Language Models (ASPO) — ICML 2025

**引用**: [Steering Protein Language Models (OpenReview)](https://openreview.net/forum?id=pu2aw2zwMx)

**核心**:
- 激活引导 (Activation Steering) 用于 PLM
- **无需训练、无需数据**
- 在热稳定/溶解度/GFP 亮度任务上有效
- 兼容 auto-encoding (ESM) 和 autoregressive (ProGen) PLM

**对本项目的应用**:
- 不重新训练,直接用 ASPO 引导 ESM2 生成
- 适用所有 GFP 候选优化

### G. InstructPLM-mu — arxiv 2025.10

**引用**: [InstructPLM-mu: 1-Hour Fine-Tuning of ESM2 Beats ESM3](https://arxiv.org/html/2510.03370v2)

**核心**:
- 1 小时 ESM2 fine-tune,性能可超过 ESM3
- 三种特征融合设计 (Cross Attention, Channel concat, Token concat)
- 适合快速集成到现有管线

**对本项目的应用**:
- ESM2-650M 1 小时微调,在 DMS 141K GFP 数据上
- 加入"亮度预测"作为新评分维度

---

## 三、GFP 设计原则 (综合 2025 文献)

### 3.1 生色团形成核心 (固定不动)
- **T65**: chromophore 形成核心
- **Y66**: chromophore 中心酪氨酸
- **G67**: 构象关键
- **R96**: 催化残基
- **E222**: 催化残基

### 3.2 折叠关键残基 (尽量保留)
- F64L (sfGFP 关键)
- F99S (折叠报告子)
- M153T (折叠关键)
- V163A (折叠关键)

### 3.3 表面优化突变 (灵活调整)
- S30R (+1.25 kcal/mol)
- Y39N, N105T, Y145F, I171V
- A206V (防二聚化)

### 3.4 设计策略层次
1. **第一层**: ESM3 提示生色团生成 (S4)
2. **第二层**: ProteinMPNN 局部优化
3. **第三层**: ThermoMPNN ΔΔG 筛选 (S3)
4. **第四层**: ESMFold 验证 + 多 recycles 投票
5. **第五层**: ESM2-3B 似然评分 (新增)

---

## 四、文献总结表

| 文献 | 发表 | 关键创新 | 对本项目价值 |
|:----|:----|:--------|:------------|
| sfGFP | NB 2006 | 11 突变 + superfolder 概念 | 基础突变库 |
| TGP | Proteins 2015 | Tm 92°C 冠军 | 极端稳定参考 |
| RFdiffusion3 | 2025.09 | 全原子扩散 | de novo 骨架 |
| Parametric β-barrel | PNAS 2025.09 | 参数化骨架生成 | 大规模 β-barrel |
| ESM3 | Science 2025 | 多模态 PLM GFP | 5 亿年进化路径 |
| ThermoMPNN | PNAS 2024 | ΔΔG 预测 | Tm 评分维度 |
| Arcadia GFP | 2025.09 | CNN ensemble | 完整 pipeline |
| ASPO | ICML 2025 | 激活引导 | 无训练优化 |
| InstructPLM-mu | arxiv 2025 | 1h 微调 | 评分维度增强 |

---

*详见 `05_next_steps.md` 的实施细节,以及 `round19/round19_next_steps.md` 的优先级*
