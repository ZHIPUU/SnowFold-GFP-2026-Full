# 📥 Round 5 需订阅账号下载的关键文献清单

> **请用您的订阅账号下载这些 PDF**，我已尽力从公开渠道提取了摘要和方法，但完整方法细节、SI Appendix、所有数据表需要全文 PDF。

---

## ⭐⭐⭐⭐⭐ 极高优先级（强烈建议下载）

### 1. EVOLVEpro — Few-shot Directed Evolution
- **期刊**: Science (Nov 21, 2024)
- **DOI**: 10.1126/science.adr6006
- **URL**: https://www.science.org/doi/10.1126/science.adr6006
- **作者**: Jiang K, Yan Z, ... Abudayyeh OO (MIT/Mass General)
- **为什么重要**: 这是 2024 年最热的"用少量实验数据 fine-tune PLM"的方法。我们的比赛**有官方 141K 数据集**，可以直接 fine-tune ESM-2，预期能极大提升突变效应预测。
- **要的内容**: 
  - Methods 完整方法
  - SI Appendix（超参数、模型选择、Active Learning 策略）
  - 各任务的具体超参表
- **bioRxiv 免费替代（早期版）**: https://www.biorxiv.org/content/10.1101/2024.07.17.604015v1

### 2. LigandMPNN — 配体感知序列设计
- **期刊**: Nature Methods (Feb 2025)
- **DOI**: 10.1038/s41592-025-02626-1
- **URL**: https://www.nature.com/articles/s41592-025-02626-1
- **作者**: Dauparas J, Lee GR, ... Baker D (UW)
- **为什么重要**: ProteinMPNN 的升级版，明确建模 chromophore (HBI/CYG/TYG) 作为配体。我们的 Round 4 用的是 ProteinMPNN 不能感知 chromophore，LigandMPNN 应该能给出更好的设计。
- **要的内容**:
  - 配体表示方法（如何把 GFP chromophore 编码）
  - 100+ 实验验证案例细节
  - SI 中的 GFP 相关案例
- **预印本免费**: https://www.biorxiv.org/content/10.1101/2023.12.22.573103v1（早期版）

### 3. ESM3 — 5 亿年进化 / esmGFP
- **期刊**: Science (Jan 16, 2025)
- **DOI**: 10.1126/science.ads0018
- **URL**: https://www.science.org/doi/10.1126/science.ads0018
- **作者**: Hayes T, Rao R, ... Rives A (EvolutionaryScale)
- **为什么重要**: esmGFP 的设计 protocol 是我们 de novo GFP 设计的金标准。58% 同源性的新荧光蛋白！
- **要的内容**:
  - GFP 设计的 chain-of-thought prompting 完整 protocol
  - 88 个候选的实验验证统计
  - SI 中 ESM3 GFP 设计的所有 prompt 细节
- **Forge API（免费试用）**: https://forge.evolutionaryscale.ai

### 4. GeoEvoBuilder — Zero-shot 稳定+活性双优化
- **期刊**: PNAS (Oct 10, 2025)
- **DOI**: 10.1073/pnas.2504117122
- **URL**: https://www.pnas.org/doi/10.1073/pnas.2504117122
- **作者**: Liu J, You H, ... Lai LH (北大来鲁华)
- **为什么重要**: GitHub 已开源代码 (`PKUliujl/GeoEvoBuilder`)，但 **SI Appendix (21.96 MB PDF)** 才有 1GFL backbone 设计的完整 hyperparameter 和 protocol。
- **要的内容**:
  - **SI Appendix.pdf**（关键！设计细节都在 SI）
  - 1GFL-15 / 1GFL-19 设计的精确突变列表
  - 晶体结构 9JVR/9JV7 的 PDB 分析

---

## ⭐⭐⭐⭐ 高优先级

### 5. Molecular Spies in Action — FP biosensor 综述（124页）
- **期刊**: Chemical Reviews (Nov 13, 2024)
- **DOI**: 10.1021/acs.chemrev.4c00293
- **URL**: https://pubs.acs.org/doi/10.1021/acs.chemrev.4c00293
- **作者**: Gest AMM, Sahan AZ, ... Zhang J (UCSD)
- **为什么重要**: FP chromophore 环境设计原则的最权威综述。124 页全是干货。
- **要的内容**: 
  - chromophore 环境 (Gln69, Arg96, Glu222 等关键残基) 的全面分析
  - 各 FP 突变的功能影响表
  - 第 2-3 章 GFP 工程原则

### 6. ProCeSa — Contrast-Enhanced Structure-Aware Thermostability
- **期刊**: J Chem Inf Model (Feb 23, 2025)
- **DOI**: 10.1021/acs.jcim.4c01752
- **URL**: https://pubs.acs.org/doi/full/10.1021/acs.jcim.4c01752
- **作者**: Zhou F, Zhang S, ... Liu JK (Birmingham)
- **为什么重要**: 基于 PLM embeddings 的稳定性预测，**不需要原子坐标**，对 GFP 很合适。
- **要的内容**: 完整方法 + 训练数据集

### 7. AI.zymes — 整合进化框架
- **期刊**: Angewandte Chemie International Edition (June 23, 2025)
- **DOI**: 10.1002/anie.202507031
- **URL**: https://onlinelibrary.wiley.com/doi/full/10.1002/anie.202507031
- **作者**: Merlicek LP, Neumann J, ... Bunzel HA (ETH Zurich)
- **为什么重要**: 整合 Rosetta + ESMFold + ProteinMPNN + FieldTools 的进化框架。
- **要的内容**: GitHub 已开源 (`bunzela/AIzymes`)，但**论文 Methods 才说明如何 tune 各模块权重**

---

## ⭐⭐⭐ 中等优先级（可选下载）

### 8. NSGA-II + ProteinMPNN 多目标优化
- **期刊**: Science and Technology of Advanced Materials: Methods (Feb 2, 2026)
- **DOI**: 10.1080/27660400.2025.2611575
- **URL**: https://www.tandfonline.com/doi/full/10.1080/27660400.2025.2611575
- **作者**: Akiba R, Moriwaki Y, ... Yoshikawa N (Tokyo)
- **为什么重要**: 多目标 ProteinMPNN, 平衡序列差异 + 结构相似度。
- **要的内容**: NSGA-II 参数 + 实现代码

### 9. AlphaDE — MCTS + PLM 定向进化
- **arXiv** (开放): https://arxiv.org/abs/2511.09900
- **作者**: Yang Y, Wang Y, ... Heng PA (CUHK / 杭州医学院)
- **为什么重要**: NeurIPS 2025 论文, 专门有 GFP benchmark
- **要的内容**: 已经在 arXiv 公开, 无需订阅 ✅

### 10. Boltz-1 — 开源 AlphaFold3
- **bioRxiv** (开放): https://www.biorxiv.org/content/10.1101/2024.11.19.624167v4
- **GitHub**: https://github.com/jwohlwend/boltz
- **为什么重要**: 替代 ESMFold 验证候选
- **要的内容**: 已开放 ✅, 直接下载 PDF

---

## 📋 已经从开放渠道获得完整信息的（无需下载）

以下已经从公开摘要、博客、GitHub README 中获得足够信息：

- ✅ RFdiffusion3 (bioRxiv + GitHub)
- ✅ SaProt (OpenReview ICLR 2024)
- ✅ ThermoMPNN (PNAS 已开放 free access)
- ✅ DCI_asym GNN (PNAS 已开放)
- ✅ htFuncLib (Nat Commun 已开放, PMC10199939)
- ✅ mBaoJin (Nat Methods 已开放, PMC11852770)
- ✅ mGreenLantern (PNAS 已开放, PMC7720163)
- ✅ Arcadia GFP CNN (Stacks 完全开放)
- ✅ ESM-IF1 (GitHub 完整文档)
- ✅ EVOLVEpro 早期版 (bioRxiv 免费)

---

## 🎯 下载建议的优先顺序

如果时间有限，按这个顺序下载就够了：

| # | 文献 | 估计阅读时间 | 立即可用价值 |
|---|------|------------|-------------|
| 1 | **EVOLVEpro Science** | 1h | ⭐⭐⭐⭐⭐ 直接用比赛数据 fine-tune |
| 2 | **LigandMPNN Nat Methods** | 1h | ⭐⭐⭐⭐⭐ 替代 Round 4 ProteinMPNN |
| 3 | **GeoEvoBuilder SI Appendix** | 2h | ⭐⭐⭐⭐⭐ 1GFL 设计 protocol |
| 4 | **ESM3 Science + SI** | 1.5h | ⭐⭐⭐⭐ esmGFP prompt 细节 |
| 5 | **Molecular Spies Chem Rev** | 跳读 4h | ⭐⭐⭐⭐ FP 原理参考书 |

---

## 💡 提取信息的备选方案

如果暂时无法下载某篇，可以尝试：

1. **bioRxiv 预印本**（许多 Science/Nature 论文有早期版）
2. **GitHub README + Wiki**（代码作者通常会写得很详细）
3. **作者实验室主页**（如 Baker lab IPD 经常发 PDF）
4. **DeepWiki**（如 ESM3: https://deepwiki.com/evolutionaryscale/esm）
5. **CSDN/知乎 中文解读**（特别是 ProteinMPNN/AlphaFold 系列）
6. **联系作者**（往往会回复 PDF 请求）

---

*请按重要度顺序下载，下载完成后告知，我将整合到 Round 5 设计中*
