# Round 5 战略规划

> **状态**: 文献调研完成，等待工具部署
> **目标**: 在 Round 4 v5 (Best Top-1 预测 1.23) 基础上，引入 2024-2026 SOTA 工具，进一步提升综合分

---

## 🎯 Round 5 核心战略

### 与 Round 4 的差异

| 维度 | Round 4 v5 | Round 5 目标 |
|------|-----------|-------------|
| **设计工具** | ProteinMPNN (vanilla) | + **LigandMPNN** (chromophore-aware) |
| **结构预测** | ESMFold (FP32 GPU) | + **Boltz-1** / ColabFold 交叉验证 |
| **突变效应预测** | 无 (放弃 ML) | + **EVOLVEpro fine-tune** (用比赛 141K 数据!) |
| **稳定性预测** | 文献先验 | + **ThermoMPNN-D** (单点 ΔΔG) + **ProCeSa** (PLM-based) |
| **零样本评分** | ESMFold pLDDT only | + **SaProt** (Structure-aware PLM) + **ESM-IF1** 评分 |
| **多目标优化** | 手动权重 | **MOG-DFM / NSGA-II** Pareto 优化 |
| **OOD 处理** | 文献+保守 | **DCI_asym GNN** 给双点突变估值 |
| **De novo 路径** | ProteinMPNN sfGFP + avGFP | + **RFdiffusion3** all-atom + **ESM3 esmGFP** + **GeoEvoBuilder** |

---

## 📋 Round 5 任务清单（按优先级）

### 🔥 P0 必做（高 ROI 即用）

#### 1. **EVOLVEpro 风格 fine-tune 比赛 141K 数据**
- **依据**: Science 2024, 用 PLM top-layer regression + Active Learning
- **数据**: 比赛官方 GFP_data.xlsx (141,572 条 CFPS 亮度)
- **流程**:
  1. 提取每条序列的 ESM-2 650M 嵌入 (Round 2 已有 `esm650m_embeddings.npy`)
  2. 训练 Random Forest top-layer regressor (论文最优)
  3. 用 trained regressor 给 53 条 Round 4 候选打分
  4. 验证：与 Round 2 XGBoost R²=0.916 对比
- **预期**: 解决 OOD 问题，给候选第二维度独立打分
- **GitHub**: https://github.com/mat10d/EvolvePro

#### 2. **LigandMPNN 替代 ProteinMPNN**
- **依据**: Nature Methods 2025
- **关键**: chromophore HBI/CYG/TYG 作为配体进入图结构
- **流程**:
  1. 准备 sfGFP/avGFP/mBaoJin PDB + chromophore 配体注释
  2. 跑 LigandMPNN 生成 100+ 候选 (与 Round 4 ProteinMPNN 对比)
  3. ESMFold pLDDT 验证
- **GitHub**: https://github.com/dauparas/LigandMPNN

#### 3. **AlphaDE GFP MCTS 复刻**
- **依据**: arXiv 2511.09900, 专门针对 GFP 设计
- **流程**:
  1. fine-tune ESM-2 在 GFP 同源序列 (用 hmmer 找)
  2. fix chromophore + β-barrel, mask 一半其他位置
  3. MCTS 搜索 + PLM fitness 引导
- **预期**: 比 ProteinMPNN 更有方向性的设计

---

### 🥈 P1 重要

#### 4. **SaProt zero-shot 突变效应评分**
- **依据**: ICLR 2024 Spotlight
- **流程**:
  1. 下载 SaProt 650M 预训练模型
  2. 用 Foldseek 提取 sfGFP/avGFP 结构 token
  3. 对每条 Round 4 候选给 zero-shot likelihood 评分
- **预期**: 比 ESM-2 更稳健的 OOD 评分

#### 5. **ThermoMPNN 单点 ΔΔG 累加预测 Tm**
- **依据**: PNAS 2024
- **流程**:
  1. 用 ThermoMPNN 跑 sfGFP/avGFP 全位点 ΔΔG
  2. 累加候选中每个突变的 ΔΔG
  3. 转换为 ΔTm (经验关系 ΔTm ≈ ΔΔG / 0.6 kcal/mol)
- **预期**: 替代我们手工 80/85/92 估值

#### 6. **Boltz-1 验证 mBaoJin 候选**
- **依据**: bioRxiv 2024.11.19, AF3 水平开源
- **流程**:
  1. 跑 Boltz-1 重新预测 M7_mBaoJin_K173R_Y196F
  2. 比对 ESMFold pLDDT 38.9 vs Boltz pLDDT
- **预期**: 验证 mBaoJin 是否真的折叠差或 ESMFold 偏见

#### 7. **htFuncLib 16K 设计数据下载**
- **依据**: Nat Commun 2023, Tm 高达 96°C
- **资源**: Addgene plasmids #191926 等系列 / 论文 SI 序列表
- **流程**:
  1. 下载 16K 设计的全序列 + 实测 Tm + 亮度
  2. 用 EVOLVEpro 框架 fine-tune
  3. 给候选打 htFuncLib 评分

---

### 🥉 P2 探索

#### 8. **RFdiffusion3 all-atom de novo**
- **关键**: 2025-12-03 release, 10× 速度
- **流程**:
  1. 输入 chromophore 配体 + 期望 β-barrel 结构
  2. 生成 100+ all-atom backbone
  3. LigandMPNN 设计序列
- **风险**: 学习曲线陡，可能成为 Round 5 主线

#### 9. **ESM3 esmGFP API 调用**
- **关键**: 申请 Forge API token
- **流程**: 复刻 Hayes 2024 的 chain-of-thought GFP 生成

#### 10. **MOG-DFM 多目标 Pareto**
- **关键**: Duke 2025-05, 离散空间直接 Pareto
- **流程**: 同时优化亮度/稳定性/折叠概率

---

## 🛠 工具实施顺序

```
Week 1 (高 ROI):
├── EVOLVEpro fine-tune 比赛 141K 数据 (P0-1)
├── LigandMPNN 跑 sfGFP + avGFP (P0-2)
└── htFuncLib 16K 数据下载 (P1-7)

Week 2 (深化):
├── SaProt zero-shot 评分 (P1-4)
├── ThermoMPNN ΔΔG 累加 (P1-5)
├── Boltz-1 验证 mBaoJin (P1-6)
└── AlphaDE GFP MCTS (P0-3)

Week 3 (探索):
├── RFdiffusion3 (P2-8)
├── ESM3 API (P2-9)
└── MOG-DFM Pareto (P2-10)

Week 4 (整合):
├── 多工具综合打分
├── Top-6 v5 重选
└── 最终提交
```

---

## 📊 预期 Round 5 收益

| 来源 | 预期 Top-1 提升 |
|------|---------------|
| EVOLVEpro fine-tune 替代猜测 | +0.2~0.3 综合分 |
| LigandMPNN 替代 ProteinMPNN | +0.1~0.2 |
| AlphaDE MCTS 优化 | +0.2 |
| SaProt 第二维度评分 | +0.1（更稳健的选择） |
| ThermoMPNN 精准 Tm | +0.1 |
| Boltz-1 拯救 mBaoJin | +0.1（如果证明 mBaoJin 实际可折叠） |

**Round 5 中性 Best Top-1 目标**: 1.5+（v4: 1.23）
**Round 5 乐观 Best Top-1 目标**: 2.5+（v4: 1.90）

---

## 📁 Round 5 目录结构（已建立）

```
work/round5/
├── literature/
│   ├── round5_literature_kb.md         # 完整文献库 (本次)
│   └── PAPERS_TO_DOWNLOAD.md           # 受订阅限制清单 (本次)
├── (待建立)
│   ├── 01_evolvepro_finetune.py        # P0-1
│   ├── 02_ligandmpnn_design.py         # P0-2
│   ├── 03_alphade_mcts.py              # P0-3
│   ├── 04_saprot_zero_shot.py          # P1-4
│   ├── 05_thermompnn_ddg.py            # P1-5
│   ├── 06_boltz_validate.py            # P1-6
│   ├── 07_htfunclib_download.py        # P1-7
│   ├── 08_rfdiffusion3.py              # P2-8
│   ├── 09_esm3_api.py                  # P2-9
│   └── 10_mogdfm_pareto.py             # P2-10
└── docs/round5/  # 由其他 Agent 撰写
```

---

## 🎯 接下来等待

1. **用户下载受订阅限制的论文**（按 PAPERS_TO_DOWNLOAD.md 顺序）
2. **用户确认 Round 5 方向**（P0/P1/P2 优先级是否调整）
3. **开始 P0 任务实施**

---

*Round 5 文献调研阶段完成 — 2026-06-22*
