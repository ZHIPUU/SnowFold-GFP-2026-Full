# 05 · 下一步方向 (Next Steps)

> **最新更新**: 2026-06-28 (R19 完成)
> **当前真实最佳**: R19 sort_score 0.9321
> **状态**: 项目接近收敛,需要换方法突破评分天花板

---

## 🎯 立刻做的事 (P0 - 提交前)

### S1. 提交 R19 (强烈推荐,5 min)
- 文件: `D:\workspace\round19\submission_r19.csv`
- 真实 sort_score = **0.9321** (项目历史最高)
- 已通过合规检查
- **6/6 全部 M 开头,239aa,不在 Exclusion_List**

### S2. R20 收敛验证 (10 min, P0 备选)
- 用 R19 Top 6 作父代,再做一轮 MPNN
- 验证是否还有提升空间
- 期望: 0.93-0.94 (边际收益递减)

---

## 🔥 短期突破 (P1 - 评分天花板)

**问题**: 生色团 pLDDT 0.95 已接近 ESMFold 评分天花板,继续 MPNN 优化无法突破

### S3. ThermoMPNN ΔΔG 预测 (1-2 h, 强烈推荐)
**来源**: [ThermoMPNN 论文 (PNAS 2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10402116/)
**原理**: 在 ProteinMPNN 基础上迁移学习,用 Megascale 数据集(180 万测量)训练
**预期**:
- 填补 Tm 预测空白
- 直接对 R19 Top 6 评估 ΔΔG
- 选择 ΔΔG < 0 (稳定化) 的候选作为最终提交
**实施**:
```bash
pip install git+https://github.com/Kuhlman-Lab/ThermoMPNN.git
# 加载预训练权重,直接推理
```
**ROI**: 中-高 (Tm 是综合分的关键维度)

### S4. ESM3 提示 GFP 生成 (4-6 h)
**来源**: [Simulating 500M years of evolution (Science 2025)](https://www.science.org/doi/10.1126/science.ads0018)
**原理**: EvolutionaryScale ESM3 (7B 参数) 可直接以生色团 6 关键残基为提示生成 GFP
**预期**:
- ESM3 生成的 GFP 与已知 GFP 序列同一性仅 58%
- "等价于模拟 5 亿年进化"
- 可能是项目最大的突破机会
**实施**:
```python
# EvolutionaryScale API (付费) 或开源 1.4B 模型
from esm.models.esm3 import ESM3
model = ESM3.from_pretrained("esm3-medium-2024-08")
```
**ROI**: 极高 (Science 论文验证)

### S5. RFdiffusion3 (RFD3) de novo β-barrel (1-2 天)
**来源**: [RFD3 论文 (2025.09)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12458353/)
**原理**: 全原子扩散模型,生成全新骨架
**预期**:
- 完全跳出 sfGFP 骨架限制
- 设计 10-20 个全新 β-barrel
- 用 ProteinMPNN 重新设计序列
**实施**:
```bash
git clone https://github.com/RosettaCommons/foundry
# 下载 RFD3 权重,运行推理
```
**ROI**: 高 (突破生色团 pLDDT 0.95 天花板)

---

## 🌱 中期探索 (P2 - 新方法)

### S6. ASPO 激活引导 (无训练) (2-3 h)
**来源**: [Steering PLMs (ICML 2025)](https://openreview.net/forum?id=pu2aw2zwMx)
**原理**: 不重新训练,直接编辑 PLM 激活引导生成
**优势**: 无需训练数据,无需微调
**实施**:
```python
from steering_plms import ASPO
aspo = ASPO(target_property="fluorescence")
steered_seq = aspo.generate(seed_seq=sfGFP_WT)
```
**ROI**: 中 (新方法,风险中等)

### S7. InstructPLM-mu ESM2 1 小时微调 (3-4 h)
**来源**: [InstructPLM-mu (arxiv 2510.03370)](https://arxiv.org/html/2510.03370v2)
**原理**: 1 小时 ESM2 fine-tune,性能可超过 ESM3
**实施**:
- ESM2-650M 在 DMS 141K GFP 数据上微调
- 用 LoRA 节省显存
- 加入"亮度预测"作为新评分维度
**ROI**: 高 (直接增强评分体系)

### S8. Arcadia GFP CNN 集成 (1-2 天)
**来源**: [Arcadia GFP 设计 (2025.09)](https://thestacks.org/publications/result-gfp-variant-design-nn)
**原理**: ESM-2 (15B) embedding + 简单 CNN ensemble 预测荧光
**价值**:
- 完整开源代码 [Arcadia GitHub](https://github.com/Arcadia-Science/2025-GFP-variant-design)
- 完整 DMS 数据集 (Zenodo)
- 实验验证流程 (RFP-GFP 双报告系统)
**ROI**: 高 (论文+代码+数据全开源)

---

## 🚀 长期 (P3 - 大型基础设施)

### S9. AlphaFold2 / ColabFold 交叉验证
- 与 ESMFold 对比,获得 ensemble pLDDT
- A800 上跑 ColabFold MSA (耗时较长)

### S10. 多模型共识评分
- 综合 ESM2 似然 + ThermoMPNN ΔΔG + ESMFold pLDDT + ESM3 likelihood
- 抗过拟合,稳健性提升

### S11. 进化搜索 + 论文突变组合
- 用遗传算法在"论文突变 + 数据驱动"组合空间搜索
- 起点:R19 Top 6,变异:每代随机选 1-2 突变替换

---

## 📊 决策矩阵

| 选项 | 估计耗时 | 期望 sort_score 提升 | 风险 | 推荐 |
|---|---|---|---|---|
| **S1 提交 R19** | 5 min | 0 (避免错误) | 极低 | ⭐⭐⭐⭐⭐ |
| S2 R20 收敛 | 10 min | +0.0-0.5% | 极低 | ⭐⭐⭐ |
| **S3 ThermoMPNN** | 1-2 h | +1-2% (Tm 维度) | 低 | ⭐⭐⭐⭐ |
| **S4 ESM3 生成** | 4-6 h | +3-5% (新骨架) | 中 | ⭐⭐⭐⭐⭐ |
| S5 RFdiffusion3 | 1-2 天 | +5-10% (突破天花板) | 中-高 | ⭐⭐⭐⭐ |
| S6 ASPO | 2-3 h | +1-3% (新方法) | 中 | ⭐⭐⭐ |
| S7 ESM2 微调 | 3-4 h | +1-2% (评分维度) | 低 | ⭐⭐⭐⭐ |
| S8 Arcadia CNN | 1-2 天 | +3-5% (完整 pipeline) | 中 | ⭐⭐⭐⭐ |

---

## 🎯 推荐策略 (按时间预算)

### 如果竞赛截止 < 1 天
1. **S1**: 直接提交 R19 ✅
2. S3: ThermoMPNN 跑一下 (1-2 h)

### 如果还有 1 周
1. S1 + S3 (Tm 维度)
2. S4 (ESM3 提示 GFP 生成)
3. S2 (R20 收敛)

### 如果还有 2 周+
1. S1 + S3 + S4
2. S5 (RFdiffusion3 de novo)
3. S7 (ESM2 微调评分)

### 最大努力方案 (2 周+)
1. **S1 + S3 + S4 + S5 + S8**
2. 期望 sort_score 突破 0.95
3. 完整对接 Arcadia 实验流程

---

## 📚 关键文献(2025 新)

1. **Parametric β-barrel design (PNAS 2025.09)** - [DOI 10.1073/pnas.2425459122](https://pnas.org/doi/10.1073/pnas.2425459122)
   - 用 RFdiffusion + 参数化生成 β-barrel
   - 高 in silico 成功率,X-ray 晶体验证

2. **RFdiffusion3 全原子 (2025.09)** - [PMC 12458353](https://pmc.ncbi.nlm.nih.gov/articles/PMC12458353/)
   - Baker Lab 2025 最新版本
   - 包含配体/核酸条件生成

3. **ESM3 5 亿年进化 (Science 2025)** - [DOI 10.1126/science.ads0018](https://www.science.org/doi/10.1126/science.ads0018)
   - 提示 GFP 生成,58% 序列同一性
   - 最直接的 GFP 设计突破

4. **ThermoMPNN ΔΔG (PNAS 2024)** - [PMC 10402116](https://pmc.ncbi.nlm.nih.gov/articles/PMC10402116/)
   - Megascale 数据集 180 万测量
   - 直接预测 ΔΔG

5. **Arcadia GFP (2025.09)** - [Arcadia GFP Design](https://thestacks.org/publications/result-gfp-variant-design-nn)
   - 完整 CNN ensemble pipeline
   - 4 个月从概念到实验验证

6. **Steering PLMs (ICML 2025)** - [OpenReview pu2aw2zwMx](https://openreview.net/forum?id=pu2aw2zwMx)
   - 无需训练的 PLM 激活引导
   - 适用于热稳定/溶解度/亮度任务

7. **InstructPLM-mu (arxiv 2510)** - [arxiv 2510.03370](https://arxiv.org/html/2510.03370v2)
   - 1 小时 ESM2 fine-tune 超过 ESM3
   - 适合快速集成到现有管线

---

## 🔄 如果截止临近 (最终方案)

**最小可行方案 (30 min)**:
1. ✅ S1: 直接提交 R19
2. ✅ S3: ThermoMPNN 跑 R19 Top 6,补充 Tm 评估

**最大努力方案 (1 周)**:
1. ✅ S1, S3
2. ✅ S4: ESM3 生成 100 条 GFP,MPNN 优化
3. ✅ S7: ESM2 1 小时微调,加入"亮度预测"维度
4. ✅ 综合排序 → 新 Top 6

---

*详见 `round19/round19_next_steps.md` 的具体实施细节,以及 `docs/06_paper_kb.md` 的论文笔记*
