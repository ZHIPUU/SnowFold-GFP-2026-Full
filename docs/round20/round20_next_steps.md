# Round 20 Next Steps

## P0 立即可做

### S1. 提交 R20 (强烈推荐, sort_score 0.9396)
- 文件: `D:\workspace\round20\submission_r20.csv`
- 真实 sort_score = **0.9396** (项目历史新高!)
- ✅ 6/6 合规通过
- **比赛截止临近**: 直接提交这个！

### S2. 等待 R22 完成 (~5 小时后)
- R22 Phase 2 (R21 MPNN) 还有 4 父代
- Phase 3 (Top 50 + r=20 重算)
- 完成后可能会更高

---

## P1 突破 0.95

### S3. ESM2-650M LoRA Fine-tune 亮度预测
**来源**: [Saadat 2025 (EPFL)](https://doi.org/10.1016/j.csbj.2025.05.022)
- 在 R11-R20 评估数据上 fine-tune
- 加入"亮度预测"作为新评分维度
- 优点: ESM2-650M 是 Hou 2025 的最优规模

### S4. ThermoMPNN ΔΔG 预测
**来源**: [ThermoMPNN (PNAS 2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10402116/)
- 等 GitHub 镜像可用时
- 填补 Tm 评分空白
- 可能让 Top 6 加分

---

## P2 长期方向

### S5. RFdiffusion3 de novo β-barrel
**来源**: [RFD3 (Baker Lab 2025.09)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12458353/)
- 完全跳出 sfGFP 骨架限制
- 用 R20 Top 6 作序列种子

### S6. ESM3 提示 GFP 生成 (Science 2025)
**来源**: [ESM3 5 亿年进化](https://www.science.org/doi/10.1126/science.ads0018)
- 以生色团 6 关键残基为提示
- 生成全新 GFP
- EvolutionaryScale API

---

## 推荐策略

**如果竞赛截止临近**:
1. **立即提交 R20** (0.9396 ✅)

**如果还有 1 周**:
1. 提交 R20 (0.9396)
2. 等 R22 完成 (可能 0.945+)
3. 部署 ESM2-650M fine-tune 进一步加分

**如果还有 2 周+**:
1. 提交 R20/R22
2. S5 (RFdiffusion3) 探索新骨架
3. S6 (ESM3) 提示 GFP
