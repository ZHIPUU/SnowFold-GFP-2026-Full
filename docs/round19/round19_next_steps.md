# Round 19 Next Steps

## P0 立即可做 (R20 - 收敛验证)

### S1. 用 R19 Top 6 再做 MPNN
- 验证是否还有提升空间(收敛判定)
- 期望: 0.93-0.94(边际收益)

### S2. ESM2-3B Fine-tune 亮度预测
- 在 A800 上跑 3B 参数 ESM2
- 微调在 141K 蛋白数据上
- 期望: 加入"亮度预测"作为新评分维度

## P1 短期 (突破评分天花板)

### S3. RFdiffusion3 de novo β-barrel
- 完全跳出 sfGFP 骨架
- 生成 10-20 个全新骨架
- 用 R18/R19 候选作为序列设计起点

### S4. ThermoMPNN ΔΔG 稳定性预测
- 填补 Tm 预测空白
- 预计能让 Tm 评分维度 +5-10°C

## P2 长期 (新方法)

### S5. 多模型共识评分
- 综合 ESM2 似然 + ThermoMPNN ΔΔG + ESMFold pLDDT
- 抗过拟合,稳健性提升

### S6. AlphaFold2 / ColabFold 交叉验证
- 与 ESMFold 对比,获得 ensemble pLDDT

---

## 推荐提交策略

**如果竞赛截止临近**:
1. **直接提交 R19** (sort_score = 0.9321, 项目最佳)

**如果还有 1 周**:
1. R20: MPNN 收敛验证 (10 分钟)
2. S2: ESM2-3B 评分 (1 天)

**如果还有 2 周+**:
1. R20 + S2 + S3 (RFdiffusion3)
2. 期望突破 0.95

---

*最后更新: 2026-06-28*
