# Round 5 文献库（Literature Knowledge Base）

> **构建日期**: 2026-06-22
> **目标**: 为 Round 5 GFP 设计提供前沿方法论支持。重点关注 2024-2026 新工具、新方法、新数据集
> **结构**: 按类别组织，每条标注 ⭐ 重要度 + 📌 是否需订阅账号下载

---

## 🎯 总体战略洞察（基于本轮全文献综述）

### 我们目前 Round 4 v5 的位置
- **MPNN_T01_014** (pLDDT 68.3, sfGFP_MPNN) 已经接近 SOTA 文献中 ESMFold "高置信度" 阈值
- **H1_avGFP_sfGFP_acid3_I152S** (Tm 92°C 估) 利用了 sf:acid.3 文献突变
- 综合分预测：1.23-1.32（中性场景）

### Round 5 可探索的新方向（按 ROI 排序）
1. **🔥 RFdiffusion3 + LigandMPNN 闭环** (Dec 2025 release) — 直接生成全新结构 + 序列
2. **🔥 EVOLVEpro few-shot active learning** (Science 2024) — 用比赛官方 14 万条数据 fine-tune ESM-2
3. **🟡 SaProt 结构感知 PLM** (ICLR 2024) — 替代 ESM-2 做 zero-shot 突变效应预测
4. **🟡 AlphaDE Monte Carlo Tree Search** (NeurIPS 2025) — fine-tune PLM + MCTS 定向进化
5. **🟢 多目标 Pareto 优化** (MOG-DFM, Duke 2025) — 同时优化亮度+稳定性
6. **🟢 Boltz-1 / AlphaFold3** — 替代 ESMFold 做更准确结构验证

---

## 📚 一类：De Novo 设计与 Inverse Folding 工具

### 1. **RFdiffusion3** ⭐⭐⭐⭐⭐ 🆕
- **发布**: 2025-12-03（IPD/David Baker lab）
- **创新**: atom-level diffusion, 10× 速度提升, 配体/DNA/小分子结合一体化
- **GitHub**: https://github.com/RosettaCommons/foundry/tree/production/models/rfd3
- **bioRxiv**: https://www.biorxiv.org/content/10.1101/2025.09.18.676967v2
- **对我们的应用**: 给定 sfGFP/avGFP backbone + chromophore 配体，生成新 sequence-structure
- **状态**: ✅ 开源, MIT 许可

### 2. **LigandMPNN** ⭐⭐⭐⭐⭐
- **Nature Methods 2025**: Dauparas et al., DOI: 10.1038/s41592-025-02626-1
- **GitHub**: https://github.com/dauparas/LigandMPNN
- **关键**: 显式建模 chromophore HBI（4-hydroxybenzylidene-imidazolinone）作为配体
- **结果**: 100+ 实验验证设计, 亲和力提升 100×
- **对我们的应用**: 比 ProteinMPNN 更适合 GFP（chromophore 作为配体让 R96/E222 等关键位更稳定）
- **PDF**: https://www.bakerlab.org/wp-content/uploads/2025/03/s41592-025-02626-1.pdf

### 3. **GeoEvoBuilder** ⭐⭐⭐⭐⭐
- **PNAS 2025-10-10**: Liu et al. DOI: 10.1073/pnas.2504117122 (北大来鲁华团队)
- **GitHub**: https://github.com/PKUliujl/GeoEvoBuilder
- **GFP 验证结果**: 设计 1GFL-15、1GFL-19 已实验确认荧光 + 热稳定提升, 晶体结构 PDB 9JVR/9JV7
- **关键**: zero-shot 同时优化稳定性和活性, 单轮可改 30% 残基
- **对我们的应用**: 直接复刻其 GFP 设计 protocol（已在 1GFL/1QY3 backbone 上验证）

### 4. **ESM-IF1** ⭐⭐⭐⭐
- **GitHub**: https://github.com/facebookresearch/esm/tree/main/examples/inverse_folding
- **关键**: 51% sequence recovery, 72% recovery for buried residues
- **对我们的应用**: 作为 ProteinMPNN 的**第二维度评分**，给候选 inverse-folding likelihood
- **示例脚本**: `score_log_likelihoods.py`

### 5. **ProteinMPNN soluble model** ⭐⭐⭐⭐
- **关键**: `--use_soluble_model` 标志,专为可溶蛋白训练
- **对我们的应用**: GFP 是可溶蛋白, 切换 soluble 模型可能给出更好的设计

---

## 📚 二类：蛋白质语言模型 (PLM) — 突变效应预测

### 6. **EVOLVEpro** ⭐⭐⭐⭐⭐ (Science 2024) 📌
- **Science 2024-11-21**: Jiang et al. DOI: 10.1126/science.adr6006
- **GitHub**: https://github.com/mat10d/EvolvePro (Abudayyeh lab)
- **创新**: PLM (ESM-2 15B) + Active Learning, 仅需 few-shot 实验数据 fine-tune
- **效果**: 5 轮即达到 160 轮 directed evolution 的效果
- **对我们的应用**: 用比赛官方提供的 **141K GFP 数据集** fine-tune ESM-2, 然后预测 6 条候选的"PLM-recalibrated"分数
- **📌 受订阅限制**: Science 论文需订阅, **建议下载**: https://www.science.org/doi/10.1126/science.adr6006
- **bioRxiv 预印本（免费）**: https://www.biorxiv.org/content/10.1101/2024.07.17.604015v1

### 7. **AlphaDE** ⭐⭐⭐⭐ (NeurIPS 2025) 🆕
- **arXiv**: https://arxiv.org/html/2511.09900v1
- **创新**: 用 MCTS + fine-tuned PLM 做定向进化, GFP 上明确测试
- **关键 GFP 实验**: 保留 chromophore + β-barrel 残基, mask 一半其他位置, 用 PLM 预测的 fluorescence fitness 引导搜索
- **对我们的应用**: 可直接套用其 GFP protocol (该论文专门设计了 GFP benchmark)

### 8. **SaProt** ⭐⭐⭐⭐ (ICLR 2024 Spotlight)
- **GitHub**: https://github.com/westlake-repl/SaProt (西湖大学)
- **创新**: 结构感知词汇 (Foldseek tokens + 残基 tokens)
- **关键**: 在 10 个下游任务全面超过 ESM-2
- **对我们的应用**: 比 ESM-2 给出更准的突变效应预测, **对 OOD 突变组合更稳健**（Round 2 失败点）

### 9. **ESM3 + esmGFP** ⭐⭐⭐⭐
- **Science 2025-01-16**: Hayes et al. DOI: 10.1126/science.ads0018
- **关键**: 生成 esmGFP, 与最近 GFP 仅 58% 同源性
- **API**: Forge https://forge.evolutionaryscale.ai
- **对我们的应用**: 用 chain-of-thought prompting + chromophore 固定生成 de novo GFP
- **deepwiki workflow**: https://deepwiki.com/evolutionaryscale/esm/6.1-gfp-design-workflow

---

## 📚 三类：稳定性/热稳定性预测

### 10. **ThermoMPNN-D** ⭐⭐⭐⭐ (PNAS 2024)
- **PNAS 2024**: Dieckhaus et al., DOI: 10.1073/pnas.2314853121
- **bioRxiv 2024.08.20.608844**: 双点突变扩展版
- **关键**: ProteinMPNN 提取嵌入 + 轻量 head, Megascale 数据集 270K 训练
- **对我们的应用**: 给每个候选的全部突变做 ΔΔG 累加（虽然双点会失败,单点ok）

### 11. **ProCeSa** ⭐⭐⭐
- **J Chem Inf Model 2025-02**: Zhou et al. DOI: 10.1021/acs.jcim.4c01752
- **关键**: PLM + 结构特征 + 对比学习, 不需要原子坐标只需 embeddings
- **对我们的应用**: 比 ThermoMPNN 更鲁棒, 可作为第二参考

### 12. **DCI_asym GNN** ⭐⭐⭐⭐ (PNAS 2025)
- **PNAS 2025**: Huynh et al. DOI: 10.1073/pnas.2502444122 (Ozkan lab, Arizona State)
- **关键**: 动力学驱动的图神经网络, 捕捉 epistasis, **专门测试 GFP brightness**
- **关键引述**: "GFP fitness corresponds to brightness, 4 distinct proteins benchmarked"
- **对我们的应用**: 解决 Round 2 的 OOD 失败问题, 给出**动力学驱动的双点突变预测**

---

## 📚 四类：结构预测工具（替代/补充 ESMFold）

### 13. **Boltz-1** ⭐⭐⭐⭐ (Dec 2024)
- **bioRxiv 2024.11.19.624167** + MIT 开源
- **GitHub**: https://github.com/jwohlwend/boltz
- **关键**: 首个完全开源达 AlphaFold3 水平
- **对我们的应用**: 替代 ESMFold 验证候选, 尤其对 mBaoJin 这种"低 pLDDT"骨架

### 14. **AlphaFold3** ⭐⭐⭐⭐
- **服务器**: https://alphafoldserver.com/ (每日有限免费配额)
- **关键**: 包含 chromophore 共因子的复合结构预测
- **对我们的应用**: 验证 chromophore 三联体是否正确形成

### 15. **ColabFold (AF2)** ⭐⭐⭐
- **GitHub**: https://github.com/sokrypton/ColabFold
- **关键**: MMseqs2 加速 MSA, 适合中等通量
- **对我们的应用**: 已成熟, 可作为 ESMFold 的二次验证

---

## 📚 五类：GFP 突变/数据集知识

### 16. **htFuncLib (sf:acid.3)** ⭐⭐⭐⭐⭐ (Nat Commun 2023)
- **Nature Communications 2023**: Weinstein et al. DOI: 10.1038/s41467-023-38099-z
- **PMC**: https://pmc.ncbi.nlm.nih.gov/articles/PMC10199939/
- **关键**: 16,000 active-site 设计, **Tm 高达 96°C**, sf:acid.3 (T65S/Q69L/S72A/T108V/Y145M/V224I) 已开源
- **Addgene plasmid**: #191926 (pET28 sf:acid.3)
- **对我们的应用**: 已用！可下载所有 16K 设计序列做 ProteinMPNN 训练补充

### 17. **mBaoJin** ⭐⭐⭐⭐⭐ (Nat Methods 2024)
- **Nature Methods 2024-04**: Zhang et al. DOI: 10.1038/s41592-024-02203-y
- **PMC**: https://pmc.ncbi.nlm.nih.gov/articles/PMC11852770/
- **8 个单体化突变**: S55T/H77R/E80G/Q140P/H141Q/C165Y/N171Y/T201A
- **晶体结构 PDB**: 8QBJ / 8Q79 / 8QDD
- **对我们的应用**: ✅ 已使用, 但 ESMFold 偏低 (39)

### 18. **mGreenLantern** ⭐⭐⭐ (PNAS 2020) 
- **PMC**: https://pmc.ncbi.nlm.nih.gov/articles/PMC7720163/
- **关键 Tm**: **87.2°C** (比 sfGFP 86.1 还高!), 抗化学/热降解
- **对我们的应用**: 可作为新 scaffold 候选, Tm > 比赛 72°C 余量 15°C

### 19. **hfYFP (hyperfolder YFP)** ⭐⭐⭐ (Nat Methods 2022)
- **PMC**: https://pmc.ncbi.nlm.nih.gov/articles/PMC9718679/
- **关键**: 抗化学/醛/锇/胍变性, 无 Cys, 氯化物不敏感
- **对我们的应用**: 它的稳定化突变可移植到 sfGFP

### 20. **Arcadia GFP CNN ensemble** ⭐⭐⭐⭐ (Sep 2025)
- **DOI**: 10.57844/arcadia-tuvv-w59k
- **GitHub**: https://github.com/Arcadia-Science/2025-GFP-variant-design
- **Zenodo 数据**: https://doi.org/10.5281/zenodo.17088256
- **关键**: CNN ensemble + 实验验证 avGFP 设计, 提供 sequence embeddings
- **对我们的应用**: 可下载预训练 CNN 直接给我们的候选打分

### 21. **CcGFP 8 (LANL)** ⭐⭐⭐ (Protein Sci 2024)
- **DOI**: 10.1002/pro.4886 (Waldo/Nguyen 团队)
- **关键**: Corynactis californica 来源, 抗 GdnHCl 极强
- **对我们的应用**: 可探索作为另一骨架

### 22. **GFP 抗体复合物极端稳定** ⭐⭐ (BMC Res Notes 2023)
- **PMC**: https://pmc.ncbi.nlm.nih.gov/articles/PMC10283196/
- **关键数据**: sfGFP Tm = **86.1°C** (CD denaturation, pH 7.0) — **重要校准!**
- **对我们的应用**: 修正了之前 78°C 的低估，72°C 加热对 sfGFP 实际余量 14°C

---

## 📚 六类：多目标优化 / 主动学习 / 整合平台

### 23. **MOG-DFM** ⭐⭐⭐⭐ (Duke 2025-05)
- **arXiv 2505.07086v2**: Multi-Objective-Guided Discrete Flow Matching
- **关键**: 在离散空间直接做 Pareto 多目标优化（不需要嵌入连续空间）
- **对我们的应用**: 同时优化亮度+稳定性+折叠概率（解决 Round 4 score 函数权重难题）

### 24. **AI.zymes** ⭐⭐⭐ (Angew Chem 2025)
- **DOI**: 10.1002/anie.202507031
- **GitHub**: https://github.com/bunzela/AIzymes
- **关键**: 整合 Rosetta + ESMFold + ProteinMPNN + FieldTools 进化框架
- **对我们的应用**: 现成的多工具集成 pipeline, 减少自己搭建工作

### 25. **NSGA-II + ProteinMPNN** ⭐⭐⭐ (Tandfonline 2026)
- **DOI**: 10.1080/27660400.2025.2611575
- **关键**: 用进化算法 + ProteinMPNN 做 Pareto 优化
- **对我们的应用**: 直接套用方法生成 100+ 序列, 然后多目标排序

---

## 📚 七类：综述类（深入学习用）

### 26. **Molecular Spies in Action** ⭐⭐⭐⭐ (Chem Rev 2024-11) 📌
- **Chemical Reviews 2024**: Gest et al. DOI: 10.1021/acs.chemrev.4c00293
- **关键**: 124 页全面回顾 FP biosensor 设计原则, chromophore 环境
- **📌 受订阅限制**: ACS Chem Rev 需订阅
- **链接**: https://pubs.acs.org/doi/10.1021/acs.chemrev.4c00293

### 27. **Biofluorsci GFP 突变综述** ⭐⭐⭐ (Jan 2026)
- **URL**: https://biofluorsci.com/posts/brighter-longerlasting-how-gfp-mutations-enhance-fluorescence-and-photostability-for-advanced-research
- **关键**: 完整 GFP 突变效应表, EGFP 突变, photostability 综述

### 28. **AI-Driven Enzyme Engineering Review** ⭐⭐⭐ (Molecules 2026)
- **PDF**: https://pmc.ncbi.nlm.nih.gov/articles/PMC12786422/
- **关键**: 整体回顾 AI 蛋白工程, AlphaFold2/RoseTTAFold/ProGen/ESM-2

### 29. **De Novo Designed Protein Toolkit** ⭐⭐⭐ (Biosafety & Health 2025-09)
- **DOI**: 10.1016/j.bsheal.2025.09.004
- **关键**: 闭环验证 + 多组学画像

---

## 📚 八类：sfGFP / CFPS 体系相关 — 重要先验

### 30. **sfGFP P. pastoris CFPS** ⭐⭐⭐
- **Frontiers Bioeng 2020**: DOI 10.3389/fbioe.2020.00536
- **关键数据**: sfGFP CFPS 产量 50.16 ± 7.49 μg/ml in 5 h batch
- **maturation rate k = 26.62×10⁻³ s⁻¹, t½ = 26.04 min**

### 31. **FAST split-GFP CFPS detection** ⭐⭐⭐ (Sci Rep 2024)
- **PMC**: https://pmc.ncbi.nlm.nih.gov/articles/PMC10997616/
- **关键 sfGFP 性质**: 半衰期 85°C ~197 s, 90°C ~73 s; refolding rate 5.0×10⁻¹ s⁻¹
- **对我们的应用**: 比赛 72°C 加热下, sfGFP 半衰期理论上**很长**（85°C 都有 3 min）

---

## 🔒 受网页拦截 / 需订阅账号下载的关键文献（按重要度排序）

> 这些是关键文献，但论文全文受订阅限制。**请用您的订阅账号下载完整 PDF**

| # | 标题 | 期刊/DOI | 重要度 | URL |
|---|------|---------|-------|-----|
| 1 | **EVOLVEpro** Few-shot directed evolution | Science 2024<br>10.1126/science.adr6006 | ⭐⭐⭐⭐⭐ | https://www.science.org/doi/10.1126/science.adr6006 |
| 2 | **LigandMPNN** (full Nat Methods version) | Nat Methods 2025<br>10.1038/s41592-025-02626-1 | ⭐⭐⭐⭐⭐ | https://www.nature.com/articles/s41592-025-02626-1 |
| 3 | **ESM3 / esmGFP** | Science 2025<br>10.1126/science.ads0018 | ⭐⭐⭐⭐⭐ | https://www.science.org/doi/10.1126/science.ads0018 |
| 4 | **GeoEvoBuilder** SI Appendix (含 GFP 设计细节) | PNAS 2025<br>10.1073/pnas.2504117122 | ⭐⭐⭐⭐⭐ | https://www.pnas.org/doi/10.1073/pnas.2504117122 |
| 5 | **Molecular Spies in Action** FP biosensor review | Chem Rev 2024<br>10.1021/acs.chemrev.4c00293 | ⭐⭐⭐⭐ | https://pubs.acs.org/doi/10.1021/acs.chemrev.4c00293 |
| 6 | **ProCeSa** Thermostability prediction | J Chem Inf Model 2025<br>10.1021/acs.jcim.4c01752 | ⭐⭐⭐ | https://pubs.acs.org/doi/full/10.1021/acs.jcim.4c01752 |
| 7 | **AI.zymes** Evolutionary enzyme design | Angew Chem Int Ed 2025<br>10.1002/anie.202507031 | ⭐⭐⭐ | https://onlinelibrary.wiley.com/doi/full/10.1002/anie.202507031 |
| 8 | **NSGA-II + ProteinMPNN** multi-objective | Sci Tech Adv Mater Methods 2026<br>10.1080/27660400.2025.2611575 | ⭐⭐⭐ | https://www.tandfonline.com/doi/full/10.1080/27660400.2025.2611575 |
| 9 | **MOG-DFM** | arXiv 2505.07086 (开放) | ⭐⭐⭐⭐ | https://arxiv.org/abs/2505.07086 |
| 10 | **Boltz-1** full paper | bioRxiv 2024.11.19 (开放) | ⭐⭐⭐⭐ | https://www.biorxiv.org/content/10.1101/2024.11.19.624167v4 |

### 优先下载顺序

1. **EVOLVEpro Science**（最核心，可立即应用到比赛 141K 数据集）
2. **LigandMPNN Nat Methods**（必读，比 ProteinMPNN 更适合 GFP）
3. **ESM3 Science**（esmGFP 设计 protocol）
4. **GeoEvoBuilder SI**（GFP 1GFL backbone 设计细节）
5. **MOG-DFM arXiv**（多目标优化, 已开放但 PDF 完整版更好读）

---

## 🚀 Round 5 行动建议

### P0 必做（高 ROI 工具立即应用）

1. **LigandMPNN**：克隆 GitHub repo，对 sfGFP 2B3P 做 chromophore-aware 设计
2. **EVOLVEpro 风格 fine-tune**：用比赛 141K 数据集 fine-tune ESM-2 的 top-layer regressor
3. **AlphaDE 的 GFP MCTS**：实现其论文中 GFP 设计 protocol（已有完整描述）
4. **htFuncLib 16K 数据下载**：Addgene plasmid 系列, 用于补充数据集

### P1 重要

5. **SaProt** 替代 ESM-2 给候选打 zero-shot 分数
6. **DCI_asym GNN** 给 OOD 双突变估值
7. **Boltz-1** 验证 mBaoJin 候选（替代 ESMFold 偏见）
8. **ESM-IF1 inverse folding** 给候选第二维度评分

### P2 探索

9. **RFdiffusion3** 完全 de novo GFP-fold
10. **AI.zymes** 整合平台快速跑闭环
11. **MOG-DFM** 多目标 Pareto 采样

---

## 📊 关键数据汇总（实验测定 FP Tm）

| FP | Tm (°C) | 来源 |
|----|---------|------|
| sfGFP | **86.1** | BMC Res Notes 2023 (CD denaturation, pH 7.0) |
| mGreenLantern | **87.2** | PNAS 2020 |
| mEGFP | ~78 | Cranfill 2016 |
| mNeonGreen | **68.0** | Campbell 2022 |
| eYFP | ~67 | Campbell 2022 |
| mBaoJin | **~92** | Nat Methods 2024 |
| StayGold dimer | **~92** | Nat Biotech 2022 |
| TGP | **85-90** | Close 2015 |
| htFuncLib top | **96** | Nat Commun 2023 |
| EGFP | **78** | Multiple |

**对比赛意义**: 72°C 处理对大多数 GFP 来说**远低于 Tm**，因此 Ffinal/Finit 应主要由
1. 加热时长（未知）
2. 复性效率
3. CFPS 体系的具体酶/盐成分
决定。

---

## 🛠 可立即使用的开源代码（GitHub）

| 工具 | URL | 用途 |
|------|-----|------|
| LigandMPNN | https://github.com/dauparas/LigandMPNN | 配体感知 sequence design |
| EVOLVEpro | https://github.com/mat10d/EvolvePro | Few-shot active learning |
| SaProt | https://github.com/westlake-repl/SaProt | Structure-aware PLM |
| GeoEvoBuilder | https://github.com/PKUliujl/GeoEvoBuilder | Zero-shot 同时优化稳定+活性 |
| Boltz-1 | https://github.com/jwohlwend/boltz | AF3-level open source |
| RFdiffusion3 | https://github.com/RosettaCommons/foundry | All-atom diffusion |
| AI.zymes | https://github.com/bunzela/AIzymes | 整合平台 |
| Arcadia GFP CNN | https://github.com/Arcadia-Science/2025-GFP-variant-design | avGFP brightness CNN |
| ColabFold | https://github.com/sokrypton/ColabFold | AF2 加速 |
| ESM (incl. IF1) | https://github.com/facebookresearch/esm | Multi-modal ESM 工具集 |
| ProteinMPNN | https://github.com/dauparas/ProteinMPNN | Round 4 已使用 |

---

*文献库构建完成 — 2026-06-22。下一步：Round 5 计划文档*
