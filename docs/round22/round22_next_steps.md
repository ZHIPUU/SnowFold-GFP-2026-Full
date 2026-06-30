# Round 22 Next Steps

## P0 立即可做 (比赛截止临近)

### S1. 提交 R22 (强烈推荐!)
- 文件: `D:\workspace\round22\submission_r22.csv`
- 真实 sort_score = **0.9430** 🏆 新项目纪录
- ✅ 6/6 合规通过
- **比赛截止临近**: 直接提交这个！

### S2. 等待 R23 (~3.5 小时)
- R23 在服务器 (`0e6f316839ca`) 跑:R20 Top 3 父代 + 高温度 × 150 × 3
- 预计 ~18:30 完成
- 完成后立即下载 + 对比 R22

---

## P1 突破 0.95+

### S3. ESM2-650M LoRA Fine-tune 亮度预测
**来源**: [Saadat 2025 (EPFL)](https://doi.org/10.1016/j.csbj.2025.05.022)
- 在 DMS 数据上 fine-tune
- "ESM2-650M 是突变预测最优规模" (Hou 2025)
- 等网络恢复

### S4. ThermoMPNN ΔΔG 预测  
**来源**: [ThermoMPNN (PNAS 2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10402116/)
- 填补 Tm 评分空白
- 等 GitHub 镜像可用
- **今天所有镜像失败,暂时搁置**

### S5. RFdiffusion3 de novo β-barrel
**来源**: [RFD3 (Baker Lab 2025.09)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12458353/)
- 用 R22 父代生成新 β-barrel 骨架
- 等 GitHub 镜像可用

### S6. ESM3 提示 GFP 生成 (Science 2025)
**来源**: [ESM3 5亿年进化](https://www.science.org/doi/10.1126/science.ads0018)
- 以 6 关键残基为提示
- 生成序列与已知 GFP 同一性仅 58%
- EvolutionaryScale API

---

## 推荐策略

**如果竞赛截止 < 1 天**:
1. ✅ S1: 提交 R22 (0.9430)
2. 仍然跑完 R23 (不会浪费)

**如果还有 1 周**:
1. S1: 提交 R22
2. S3: ESM2-650M Fine-tune (网络恢复后)
3. S5: RFdiffusion3 de novo 骨架
4. S6: ESM3 提示 GFP

**最大努力方案 (2 周+)**:
1. S1 + S2 + S3 + S5 + S6
2. 期望突破 0.96

---

## 中期目标

| 分数 | 突破方式 |
|:----:|:--------|
| **0.93** | ✅ R20 fixed position 修复后达成 |
| **0.94** | ✅ R22 Phase 2 大规模 MPNN 达成 |
| **0.95** | 需要 ESM2-650M Fine-tune 评分维度 |
| **0.96+** | 需要 RFdiffusion3 / ESM3 新骨架 |

---

## 当前最佳 Top 6 (R22)

```
Seq1 = 0.9430  Seq2 = 0.9429  Seq3 = 0.9416
Seq4 = 0.9416  Seq5 = 0.9414  Seq6 = 0.9413
```

文件: `D:\workspace\round22\submission_r22.csv`

---

*最后更新: 2026-06-29 14:21*
