# SnowFold GFP Design 2026 — 完整实验记录

> 从公开 PDB 结构出发，经 27 轮迭代设计，sort_score 从 0.80 提升至 **0.9477**。

## 项目简介

本仓库包含 2026 合成生物 GFP 设计竞赛的**完整实验记录**，涵盖从 R2 到 R27 的全部脚本、中间文件、MPNN 输出、分析工具和文档。

## 实验数据全景数据库与可视化平台

本项目构建了一套完整的实验数据管理系统，将 27 轮迭代的全部数据整合为 SQLite 数据库 + 交互式网站。

### 数据库规模

| 表 | 记录数 | 内容 |
|:---|------:|:-----|
| sequences | 243,386 | 去重序列（SHA-256 主键） |
| metrics | 1,474,887 | 全部预测指标 |
| rounds | 27 | 轮次元数据 |
| artifacts | 5,539 | 文件证据索引 |
| lineage_edges | 16,000+ | 序列谱系关系 |
| documents | 157 | 实验文档全文 |

### 可视化平台

- **技术栈**: FastAPI + React + Cytoscape.js + Plotly + R/ggplot2 + Docker
- **功能**: Dashboard / Sequence Vault / Topology Network / Docs Browser / R 图表报告
- **位置**: `protein-design-atlas/`

## 关联仓库

本项目由三个仓库组成：

| 仓库 | 说明 | 链接 |
|:-----|:-----|:-----|
| **竞赛核心仓库** | 标准化 pipeline 脚本 + 最终结果 | [SnowFold-GFP-2026](https://github.com/ZHIPUU/SnowFold-GFP-2026) |
| **完整实验记录** (本仓库) | R2-R27 全部脚本/数据/文档 | [SnowFold-GFP-2026-Full](https://github.com/ZHIPUU/SnowFold-GFP-2026-Full) |
| **gssh CLI** | 远程 GPU 服务器管理工具 | [gssh](https://github.com/ZHIPUU/gssh) |
| **Protein Design Atlas** | 数据库 + 交互式可视化平台 ([在线访问](http://120.48.98.164:18082)) | [protein-design-atlas](https://github.com/ZHIPUU/protein-design-atlas) |

> 由于本项目跨越 27 轮迭代实验，原始工作目录中包含数百个脚本文件（含早期 ML 预测、多版本调试、中间筛选步骤、MPNN FASTA 输出、PDB 结构等），即使按轮次分类后每轮的相关文件仍然繁多且结构不一。因此我们另建了竞赛核心仓库，仅保留标准化后的核心 pipeline 脚本和最终结果。本仓库为完整原始记录。

## 最终成绩

| 指标 | 值 |
|:-----|:---|
| **sort_score (Top-1)** | **0.9477** |
| pTM | 0.9321 |
| 全局 pLDDT | 94.65 |
| 生色团 pLDDT | 96.96 |
| 提交轮次 | R25 |

## 完整设计链条

```
公开 PDB (2B3P sfGFP / 2WUR avGFP / 8QBJ mBaoJin)
  ↓ R2: XGBoost + ESM2-650M 预测
  ↓ R4: ProteinMPNN 多骨架设计
  ↓ R5-R14: 迭代 MPNN + ESMFold
  ↓ R17-R18: 标准化 pipeline
  ↓ R19-R20: fixed pos 1+M 修复 → 0.9396
  ↓ R22: Phase 1+2 大规模 → 0.9430
  ↓ R24: 跨度温度 → 0.9447
  ↓ R25: 中低温精细 → 0.9477 🏆
```

## 目录结构

```
2026Protein Design/
├── docs/                    # 各轮实验文档 (R3-R27)
├── work/                    # 早期实验工作目录 (R2-R8, 含 ProteinMPNN/LigandMPNN)
├── scripts/                 # 标准化脚本
│   ├── pipeline/            # R17-R28 pipeline 脚本
│   └── analysis/            # 分析与校验工具
├── results/                 # 各轮 Top6 结果 (R17-R27)
├── protein-design-atlas/    # 数据库 + 可视化网站 (SQLite + FastAPI + React)
├── referencepaper/          # 参考论文 (StayGold, mBaoJin, superfolder)
├── Exclusion_List.csv       # 竞赛排除序列表
├── GFP_data.xlsx            # 竞赛提供的 GFP 突变数据
└── README.md                # 本文件
```

## 自主开发工具：gssh CLI

本项目全程使用了我们自主开发的 **gssh** 命令行工具——一个专为远程 GPU 服务器协作而设计的轻量 CLI。

- **仓库地址**: [https://github.com/ZHIPUU/gssh](https://github.com/ZHIPUU/gssh)
- **功能**: 远程任务执行、文件传输、实时日志、后台任务管理
- **应用**: R20-R27 全部远程任务的启动、监控、结果下载均通过 gssh 完成

## 技术栈

- ProteinMPNN (v_48_020) — 序列生成/反向折叠
- ESMFold (facebook/esmfold_v1) — 结构预测
- ESM2-650M — 序列嵌入
- XGBoost — 早期亮度预测
- gssh — 远程 GPU 服务器管理
- Trae AI Agent (Claude) — 辅助设计与迭代

## License

MIT License
