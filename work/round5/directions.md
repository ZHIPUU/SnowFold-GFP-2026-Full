# Round 5 剩余可探索方向（按 ROI 排序）

> **当前状态**: Round 5 v1 完成，Best Top-1 预测 1.23 (与 Round 4 v5 持平)
> **核心瓶颈**: 预测分受限于"突变数越多越不安全"的保守假设，而 LigandMPNN 的 170 突变高 pLDDT 候选没有被充分奖励

---

## 已完成 vs 未完成 总览

| 工具 | 状态 | 实际效果 |
|------|------|---------|
| ProteinMPNN (sfGFP) | ✅ Round 4 完成 | pLDDT 68.3 (MPNN_T01_014) |
| ProteinMPNN (avGFP) | ✅ Round 4 完成 | pLDDT 61.5 (MPNN_av) |
| LigandMPNN (avGFP) | ✅ Round 5 完成 | pLDDT 63.2 (R5_av_lmpnn_v2_025) |
| LigandMPNN (sfGFP) | ❌ CRO chromophore 问题 | 需修复 |
| GeoEvoBuilder | ❌ 缺权重文件 | 需北大网盘下载 |
| EVOLVEpro | ❌ 未开始 | 高 ROI |
| ThermoMPNN | ❌ 未开始 | 中 ROI |
| SaProt | ❌ 未开始 | 中 ROI |
| Boltz-1 | ❌ 未开始 | 中 ROI |
| ESM3 Forge API | ❌ 未申请 | 低 ROI |
| RFdiffusion3 | ❌ 未开始 | 探索性 |

---

## 🔥 P0：(1-2 小时完成，潜在收益最大)

### 1. 修复 sfGFP LigandMPNN (CRO → TYG)
**问题**: 2B3P 的 chromophore 是 CRO (HETATM)，LigandMPNN 把这个当成可重设计位点
**方案**: 用 afGFP 的"原生 TYG"PDB（如 1QY3 或 2WUR-parsed），或直接用我们已有的 avGFP PDB 映射到 sfGFP 序列
**预期**: 生成 30+ 条 sfGFP chromophore-retained 候选，pLDDT 可能 > 65
**时间**: 30 分钟

### 2. 用比赛 141K 数据做 EVOLVEpro 风格 fine-tune
**方案**: 
- 用 Round 2 已有的 `esm650m_embeddings.npy` (141K × 1280)
- 训练 Random Forest top-layer (EVOLVEpro 论文最优)
- 给 71 条候选做独立评分
**预期**: 解决 OOD 问题，给候选第二维度分数
**时间**: 1-2 小时

### 3. GeoEvoBuilder 直接复刻 1GFL 设计
**方案**: 论文已公开 GFP 设计 protocol，1GFL backbone 下载完成
**需要**: 北大网盘权重 `Se.pt` (访问码 `xx7W`)
**预期**: 零样本同时优化荧光+热稳定，论文实验已确认效果
**时间**: 权重下载后 30 分钟

---

## 🟡 P1：(2-4 小时，稳健提升)

### 4. ThermoMPNN 单点 ΔΔG 累加 Tm 预测
**方案**: 在 codeGPT 上找到 ThermoMPNN 推理代码，跑 sfGFP/avGFP 全位点 ΔΔG
**预期**: 替代手工 80/85/92 估值，更精准 Tm
**时间**: 2 小时

### 5. SaProt 零样本突变效应评分
**方案**: 下载 SaProt 650M 模型，用 Foldseek 提取结构 token
**预期**: 比 ESM-2 更稳健的 OOD 评分
**时间**: 3 小时

### 6. Boltz-1 验证 mBaoJin 候选
**方案**: 用 Boltz-1 替代 ESMFold 重新预测 M7_mBaoJin_K173R_Y196F
**预期**: 验证 mBaoJin 是否真的折叠差或 ESMFold 偏见
**时间**: 2 小时

### 7. 多目标优化框架 (MOG-DFM / NSGA-II)
**方案**: 对 71 条候选做 Pareto 排序，同时优化亮度+稳定性+折叠概率
**预期**: 更科学的 Top-6 选择
**时间**: 2 小时

---

## 🟢 P2：探索性（收益不确定）

### 8. ESM3 Forge API — esmGFP 设计
**方案**: 申请 Forge API token，复刻 Hayes 2024 的 chain-of-thought GFP 生成
**限制**: a) API 可能需要付费 b) 网络可能在境外
**时间**: 4 小时

### 9. RFdiffusion3 — all-atom de novo
**方案**: 下载 RFdiffusion3 模型，给定 chromophore 配体生成全新 backbone
**限制**: 可能需要 32GB+ VRAM，RTX 5080 16GB 可能不够
**时间**: 1 天

### 10. AlphaDE MCTS 定向进化
**方案**: 实现 arXiv 2511.09900 的 GFP 设计 protocol
**限制**: 需要 fine-tune ESM-2，计算量大
**时间**: 2-3 天

---

## 📝 比赛必交（必须完成）

### 11. 设计思路 PDF 文档
**要求**: 展示完整 pipeline、算法管线、6 条候选选择原因
**时间**: 2-3 小时

### 12. GitHub 开源仓库
**要求**: 完整代码、README、环境配置、一键运行
**时间**: 1 小时

---

## 🎯 推荐执行顺序

```
今天 (2-3 小时):
├── P0-1: 修复 sfGFP LigandMPNN (30 min)
├── P0-2: EVOLVEpro fine-tune 141K 数据 (1-2 h)
└── P0-3: GeoEvoBuilder (需权重下载)

明天:
├── P1-4: ThermoMPNN ΔΔG (2 h)
├── P1-7: 多目标 Pareto 优化 (2 h)
└── P0-2 续: EVOLVEpro 完整评估

本周内:
├── P1-5: SaProt 评分 (3 h)
├── P1-6: Boltz-1 验证 (2 h)
└── 必交: 设计 PDF + GitHub 仓库
```

---

## 🚨 最可能产生突破的方向（按概率）

| 方向 | 突破概率 | 原因 |
|------|---------|------|
| **sfGFP LigandMPNN 修复** | 80% | 已有 avGFP 成功经验，只需换 PDB 格式 |
| **EVOLVEpro fine-tune** | 70% | 比赛 141K 数据是现成的，ESM-2 嵌入已有 |
| **GeoEvoBuilder** | 60% | 论文实验验证过 GFP，但权重需要下载 |
| **ThermoMPNN** | 80% | 成熟工具，但 ΔΔG 累加对多突变失效 |
| **Boltz-1** | 50% | 可能给出不同 pLDDT，但不一定更好 |
| **ESM3 API** | 40% | API 限制 + 付费风险 |

---

## 💡 我的建议

**最高 ROI 的三件事**（今天就能完成）：

1. **修复 sfGFP LigandMPNN** — 30 分钟，pLDDT 可能 > 65，直接加入提交
2. **EVOLVEpro 风格 fine-tune** — 用比赛 141K 数据训练顶层回归器，给 71 候选打分，这是**唯一能解决 OOD 问题的方法**
3. **GeoEvoBuilder** — 如果你能下载北大网盘的 `Se.pt` 权重文件

之后再做 ThermoMPNN + 多目标优化，最后写 PDF + GitHub。

你觉得这个优先级如何？想从哪个方向开始？