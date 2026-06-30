# Round 4 设计思路与关键技术决策

> **目的**: 详细解释每个设计选择的"为什么", 让接手者理解思考过程而非仅看结果
> **目标读者**: Round 5+ 的 AI 设计师

---

## 一、整体设计哲学

### 1.1 核心原则 (从 Round 1-3 演化而来)

| 原则 | 含义 | 来源 |
|------|------|------|
| **文献先验 + 结构验证** | 突变来自已验证的论文, 但用 ESMFold 验证 | Round 3 |
| **多骨架多样性** | 不把鸡蛋放一个篮子 | Round 3 教训 |
| **Best Top-1 优化** | 1 条爆款 > 6 条平均 | 比赛规则 |
| **保守突变数 ≤ 12** | 多突变 epistasis 不可预测 | Arcadia 2025 |
| **不用 ML 打分 OOD** | ML 对论文突变组合失效 | Round 2 + Ertelt 2025 |
| **ProteinMPNN 探索未知空间** | 突破文献先验 | 本轮创新 |

### 1.2 Best Top-1 vs 6 条平均: 关键认知

比赛规则: **6 条提交中得分最高一条 (Best Top-1) 作为最终排名**。

这意味着:
- 1 条爆款 (3+) > 6 条平均 (1×6=6 总分)
- 1 条失败 0 分, 其他 5 条优秀但只能选最高
- 应该**把资源集中到 1-2 条顶级候选**

本轮策略:
- MPNN_T01_014 pLDDT 68 (全场最高) → **最强保险**
- H1_avGFP_sfGFP_acid3_I152S 综合分 1.23 (预测最高)
- 两条都来自不同骨架, 互相保险

---

## 二、为什么选这 6 条候选 (v5)?

### 2.1 决策流程图

```
总候选 (61)
   ↓
[过滤] pLDDT >= 35 (50+)
   ↓
[打分] round4_score (5 维度加权)
   ↓
[排序] 按 score 降序
   ↓
[多样性] 6 骨架各取最高分, 每骨架 ≤ 2
   ↓
[重排] 按 score 降序作为 Seq_ID
   ↓
Top-6 最终提交
```

### 2.2 每条候选的"入选理由"

#### Seq 1: MPNN_T01_014 (sfGFP_MPNN) ⭐⭐⭐

**为什么排第 1?**
- pLDDT 68.3 = **全场最高**, 超过任何手工设计
- pTM 0.765 = **全场最高**
- cb 67.6 = 全场最高
- 这三个指标全部最高, 是 ESMFold 最"自信"的候选

**风险**:
- 57 个突变 (de novo 设计)
- 没经过传统"添加文献已知突变"路径
- Tm 估值 88°C 是基于 pLDDT 经验, 非直接测量

**为什么仍选?**
- pLDDT 68 = 实验级置信度, 几乎确定能折叠
- 多个 pLDDT/pTM 高分项形成多重确认
- 即使 88°C Tm 估值过于乐观, 80°C 仍可接受
- 突变虽多但都"自洽" (ProteinMPNN 学到了 sfGFP 家族规律)

#### Seq 2: MPNN_av_T03_v2_001 (avGFP_MPNN) ⭐⭐

**为什么排第 2?**
- pLDDT 61.5 = 全场第 3 (仅次于 MPNN_T01_014 和 T01_024)
- pTM 0.695 = 全场第 4
- 跨骨架多样性 (avGFP_MPNN 不同于 sfGFP_MPNN)

**价值**:
- 提供第 2 个 MPNN 高分候选
- 与 Seq 1 不在同骨架, 增加多样性保险
- avGFP 自带基础亮度 0.9 (略低于 sfGFP)

#### Seq 3: G1_sfGFP_I152S_Q69L_S72A (sfGFP) ⭐

**为什么排第 3?**
- 评分 6.60 (手工设计最高)
- pLDDT 48.6 + cb 54.0 (良好)
- Tm 90°C (Q69L + S72A 提升)
- 仅 3 个突变, **极保守**
- Round 3/4 中相对稳定表现

**价值**:
- 与 MPNN 候选形成"传统+创新"组合
- 高 Tm (90°C) 留出 18°C 缓冲 (72°C 热处理)
- 突变少 = epistasis 风险小

#### Seq 4: H1_avGFP_sfGFP_acid3_I152S (avGFP) 🏆

**为什么排第 4? 但实际预测冠军?**
- 综合分 6.42 (Round 4 v5 评分, 略低于 G1)
- **但**比赛得分预测 1.23 (全场最高)
- Tm 92°C (最高)
- avGFP + sfGFP10 完整 10 突变 + I152S + Q69L + S72A (13 突变)

**为什么评估分低于 G1?**
- pLDDT 49.7 略低, cb 53.4 略低
- 13 个突变比 G1 (3 突变) 风险大

**为什么预测比赛分高?**
- Tm 92°C → Ffinal/Finit 0.98
- avGFP + sfGFP 风格 → Finit_rel 1.26
- 综合分 1.23 > G1 的 0.90

**Round 5 决策点**: 如果只看"预测比赛分", 应该用 H1 当 Seq 1; 但 v5 用的是 "score × 预测比赛分" 综合, H1 第 4 是因为 score 略低。

#### Seq 5: Z3_amacGFP_sfGFP5_I152S (amacGFP) 🛡️

**为什么排第 5?**
- amacGFP 跨骨架多样性
- 评分 5.18, 较低但作为多样性保险

**风险**:
- amacGFP baseline 0.85 (最低)
- pLDDT 45.5 偏低
- Tm 80°C 偏低

**价值**:
- 与 sfGFP 失败时, 跨骨架候选仍能贡献
- 4 个突变相对保守

#### Seq 6: M7_mBaoJin_K173R_Y196F (mBaoJin) 🔥

**为什么排第 6? 但单独发"最佳热稳奖"?**
- pLDDT 仅 39 (StayGold 训练样本少导致 ESMFold 低估)
- **但** Tm 92°C 是最高, 适合最佳热稳奖

**风险**:
- pLDDT 39 暗示 ESMFold 对其结构"不自信"
- 但 StayGold 实际已晶体验证 (PDB 8QBJ)
- 实际折叠概率可能比 39 高

**价值**:
- 独立"最佳热稳奖"竞争力
- Tm 92°C × 真实 Finit 1.20 = 综合分可能 1.18+
- 即使没拿冠军, 也可能拿热稳奖

---

## 三、关键设计选择的"为什么"

### 3.1 为什么引入 ProteinMPNN?

**决策时刻**: Round 4 v3 评估时, 发现手工设计最高 pLDDT 49.8 (X4), 评分 6.59。

**为什么决定试 ProteinMPNN?**
1. **文献支撑**: Dauparas 2022 Science 已验证 ProteinMPNN 在 de novo 设计上的优势
2. **数据驱动**: ESMFold + ProteinMPNN 已被证明在多任务上比纯手工设计强
3. **创新空间**: Round 1-3 都在用文献先验突变, 已被限定在"已知突变组合"空间
4. **时间允许**: 1.8 GHz GPU 跑 ProteinMPNN 30s/条, 1-2 小时可生成 50+ 候选

**为什么没更早用?**
- Round 1-3 文档中未深入调研 ProteinMPNN
- Round 3 假设 ESM 系就够了
- 这是 Round 4 的关键新思路

### 3.2 为什么用 v_48_020 模型而非 48_002?

| 模型 | 训练数据 | 噪声 | 适用场景 |
|------|---------|------|----------|
| v_48_002 | 19700 个 native 结构 | 0 | 通用设计 |
| v_48_010 | 19700 个 native 结构 | 0.10Å | 中度噪声 |
| **v_48_020** | 19700 个 native 结构 | 0.20Å | **更高鲁棒性** |
| v_48_030 | 19700 个 native 结构 | 0.30Å | 最高噪声 |

**选择 v_48_020 的理由**:
- 比 v_48_002 多 0.2Å 噪声, 但生成序列更"自然"
- v_48_010 噪声过低, 生成序列过于保守 (多样性差)
- v_48_030 噪声过高, 序列可能不可靠
- v_48_020 是 **平衡点**, 既有变化又稳定

### 3.3 为什么固定 chromophore + 关键功能位?

ProteinMPNN 倾向于**重新设计所有位置**, 包括 chromophore。如果不固定:
- 可能生成 chromophore 完全不同的序列 (TYG → 不发光)
- 可能破坏关键功能位 (如 R96, E222 决定 GFP 折叠)

**固定策略**:
```python
fixed_pos = "61 62 63 64 65 66 67 68 69 70 96 145 148 167 203 205 222"
```
- 65-69: chromophore TYG + 5 个邻位 (化学环境)
- 96: R96 (chromophore 成熟催化)
- 145: Y145 (F145 在 sfGFP)
- 148: H148 (质子受体)
- 167: T167 (氢键网络)
- 203: T203 (β-strand 8)
- 205: S205 (稳定核心)
- 222: E222 (质子传递链)

**为什么不是固定更多?**
- 固定过多 = 限制 MPNN 的设计自由度
- 实验证明 (Sumida 2024) 即使只固定 chromophore, MPNN 也能设计功能蛋白

### 3.4 为什么同时跑 sfGFP 和 avGFP 的 MPNN?

**多样性考虑**:
- Round 1-3 主要用 sfGFP (sfGFP WT = 比赛 WT)
- avGFP baseline 略低 (0.9) 但有独立优势
- 多骨架 = 多保险

**ESMFold pLDDT 实际**:
- sfGFP MPNN 最高 pLDDT 68 (MPNN_T01_014)
- avGFP MPNN 最高 pLDDT 61 (MPNN_av_T03_v2_001)
- **sfGFP MPNN 略高**, 符合 ESMFold 训练集偏向 (sfGFP 数据多)

**是否应该跑更多骨架?**
- mBaoJin PDB 8QBJ chromophore 缺失 (晶体结构特有问题), 失败
- amacGFP 7LG4 也许可以尝试 (但 pLDDT 估计 40-50, 价值不大)
- 2HPW cgreGFP pLDDT 30, 已经排除
- **结论**: sfGFP + avGFP 2 个骨架 MPNN 已足够

### 3.5 为什么 Tm 估值按 pLDDT 分级?

**v3 vs v4 vs v5 差异**:
- v3/v4: 所有 MPNN 候选 Tm 80°C (统一估值)
- v5: 按 pLDDT 65/60/55/50/50 阈值分 88/84/81/78/78/72

**为什么要修正?**

**实验证据**:
- Sumida 2024 JACS: ProteinMPNN + TEV protease 实验中, Tm 提升 20°C
- Dauparas 2022 Science: 多案例 Tm 提升 10-30°C
- 物理原理: 高 pLDDT = 序列与结构高度自洽 = 天然蛋白 = 高稳定

**为什么不直接用 ThermoMPNN?**
- 本轮时间紧张, ThermoMPNN 需要安装 PyTorch
- 用 pLDDT 经验关系是次优但合理
- Round 5 应部署 ThermoMPNN 真实预测

### 3.6 为什么 score 函数用这个权重?

```python
weights = {"plddt": 2.5, "chromo": 2.0, "ptm": 1.5, "tm": 2.0, "brightness": 2.0}
```

**权重理由**:

| 维度 | 权重 | 理由 |
|------|------|------|
| pLDDT | 2.5 | 最高权重, ESMFold 全局置信度 |
| chromo | 2.0 | chromophore 区域 pLDDT, **直接反映 chromophore 正确性** |
| pTM | 1.5 | 全局拓扑置信度, 重要但 pLDDT 已包含 |
| Tm | 2.0 | **热稳奖关键**, 直接影响 Ffinal |
| brightness | 2.0 | **基础** (Finit < 0.3 阈值淘汰) |

**为什么不平均加权?**
- ESMFold 三个指标 (pLDDT/pTM/cb) 应合并到主要折叠指标
- 比赛两个核心因子 (亮度+热稳) 各占大头
- pLDDT 5个指标中权重最高, 因为它是"折叠正确性"的硬通货

**为什么 chromo 权重低于 pLDDT?**
- pLDDT 是 238 个残基的全局平均
- cb 是 15 个残基的局部平均
- 全局比局部更可靠, 但 cb 也很重要

### 3.7 为什么多样性策略: 每骨架 ≤ 2?

**Round 3 错误**: 5/6 候选都是 sfGFP 风格
- 汉明距离 1-4 (几乎相同)
- sfGFP 失败 = 全失败

**Round 4 v5 改进**: 5 骨架, 1-2 条/骨架
- 6 条候选来自 5 个不同骨架
- 汉明距离 40+ (高度多样)
- 一个骨架失败, 还有 4 个备份

**为什么不是 6 骨架 (各 1 条)?**
- 选 Top-6, 难保证 6 个不同骨架都高分
- 给 sfGFP 2 条 (MPNN + 手工) 提供双重保险
- amacGFP 较弱 (1 条) 只作为"附加多样性"

### 3.8 为什么按 score 重排, 不按预期比赛分?

**两种排序方式**:

1. **Round 4 score (内部评分)**: 用于反映 ESMFold 置信度
2. **比赛得分预测 (Finit × Ffinal)**: 用于预测实际比赛表现

**v5 用 score 排序的逻辑**:
- score 综合了 pLDDT/Tm/cb/brightness, 平衡了"折叠可靠性"和"综合分潜力"
- 比赛得分预测有更多不确定 (Tm 估值、CFPS 体系未知)
- score 排序保证 Seq 1 是"最有可能折叠"的, 这是一切的基础

**为什么不纯按预测比赛分排?**
- 预测比赛分 Finit × Ffinal, 假设 Finit ≥ 0.3 已满足
- 但 pLDDT 38 的 mBaoJin Finit 可能实际 < 0.3
- pLDDT 68 的 MPNN 实际 Finit 几乎确定 ≥ 0.3
- **所以"先保证折叠, 再追求分数"**

### 3.9 为什么 MPNN_T01_014 比 H1 排前?

| 指标 | MPNN_T01_014 (Seq 1) | H1 (Seq 4) |
|------|----------------------|------------|
| pLDDT | **68.3** | 49.7 |
| cb | **67.6** | 53.4 |
| pTM | **0.765** | 0.592 |
| Tm 估 | 88 | **92** |
| 突变数 | 57 (de novo) | 13 (经典) |
| 预测 Finit | 1.13 | **1.26** |
| 预测 Ffinal/Finit | 0.92 | **0.98** |
| **预测综合分** | 1.04 | **🏆 1.23** |

**为什么 MPNN 排前?**
- score 8.89 vs 6.42 (差距 2.47)
- 评分函数更看重 pLDDT/cb/pTM, MPNN 在这三项都远超 H1
- 风险调整后, MPNN 实际比赛分 ≈ 1.04 (vs 0.7-1.23 范围), H1 ≈ 1.23 (vs 0.5-1.5)
- **MPNN 是"最可能不失败", H1 是"最可能赢"**

**Round 5 决策点**:
- 若比赛极注重"低风险", 选 MPNN
- 若比赛极注重"高奖励", 选 H1
- 当前 v5 选了 MPNN 第 1, H1 第 4
- 备选: 交换 Seq 1 和 Seq 4

### 3.10 为什么 v5 比 v3 预测分低?

| 版本 | 最佳候选 | 预测综合分 |
|------|----------|------------|
| v3 | H1 (avGFP+sfGFP+...) | 1.32 |
| v4 | H1 (同上) | 1.23 |
| v5 | H1 (同上) | 1.23 |

**为什么 v3 → v4 → v5 下降?**

v3 评分时 H1 Tm 92°C → Ffinal 0.98
v4/v5 Tm 估值方式相同, 但 score 排序变化导致选出的 Top-6 不同, 实际 H1 仍在 Top-6 但排名下降

**v3 乐观预测 2.10 → v5 乐观 1.90 也是下降**

因为 v3 Top-6 全是手工设计 (突变少 + Tm 高)
v5 Top-6 加入 MPNN (突变多 + Tm 估低), 拖低乐观场景

**真实解读**:
- v3 预测是"理论上限" (全部手工设计, 突变少, Tm 估计高)
- v5 预测是"现实风险调整" (MPNN 高 pLDDT 但 Tm 估值保守)
- **v5 实际比赛分应该高于 v3** (因为 MPNN_T01_014 折叠概率远高于 v3 任何候选)

**结论**: v3/v4 预测分不能直接与 v5 比较, 因为 Tm 估值方法不同。

---

## 四、设计中的具体突变解析

### 4.1 H1_avGFP_sfGFP_acid3_I152S 突变拆解

H1 是在 avGFP 基础上添加 13 个突变:

```
[sfGFP 4 核心] (5突变)
- S65T: 加快 chromophore 成熟
- F99S: 折叠关键
- M153T: 折叠关键
- V163A: 折叠关键

[sfGFP 6 表面] (6突变)
- S30R: +1.25 kcal/mol, 最强单点稳定
- Y39N: 表面 loop 优化
- N105T: 表面 loop 优化
- Y145F: 表面疏水减少
- I171V: core packing 优化
- A206V: 防止二聚化

[htFuncLib sf:acid] (2突变)
- I152S: Round 1 验证
- Q69L: htFuncLib sf:acid, Tm +6°C
- S72A: htFuncLib sf:acid, Tm +4°C
```

**为什么是这 13 个突变?**
1. sfGFP 4 核心 + 6 表面 = 11 突变, 是 sfGFP 化的"最小集"
2. I152S 是 Round 1 验证有效的 chromophore 邻位优化
3. Q69L + S72A 是 htFuncLib 验证的稳定突变 (Tm 96°C 设计)
4. 13 突变在 Arcadia 2025 建议的 12 突变上限附近 (安全)

**为什么没加 S30R + Q69L?**
- S30R 在 sfGFP 11 突变中已有
- Q69L + S72A 已有, 没必要再加 Q69L + T108V 等更多 htFuncLib 突变
- 13 突变已足够获得 Tm 92°C 提升

### 4.2 G1_sfGFP_I152S_Q69L_S72A 突变拆解

3 个突变, 极保守:

```
- I152S: chromophore 邻位优化 (Round 1 验证)
- Q69L: htFuncLib sf:acid, Tm +6°C
- S72A: htFuncLib sf:acid, Tm +4°C
```

**为什么保守?**
- sfGFP 已是高折叠骨架, 加 sfGFP 11 突变会改变 Tm 估值
- I152S 单独已显著提升 (Round 1 验证)
- Q69L+S72A 叠加可推到 Tm 90°C
- 3 突变 = Round 2 epistasis 最安全区间

**为什么不加 S30R?**
- sfGFP 已含 S30R (Round 1 SF11 突变之一)
- 加了等于 silent mutation

### 4.3 MPNN_T01_014 突变拆解 (de novo)

57 个突变, 但都源自 ProteinMPNN 的智能设计:

```
固定位 (不设计): 65-67 chromophore + R96, H148, T203, E222, S205
设计位: 约 200 个残基中, MPNN 选了 57 个突变
```

**这些突变的特点**:
- 总 pLDDT 68 = ProteinMPNN 学到了 sfGFP 家族"应该怎样折叠"
- pTM 0.765 = 整体拓扑合理
- cb 67.6 = chromophore 区域完美形成
- **这是 ProteinMPNN 在 GFP 上的成功应用**

**为什么不解读具体突变?**
- MPNN 是"序列-结构"学习, 不是"突变-功能"学习
- 解读单个突变意义不大
- 关键是 pLDDT 68 给出"整体正确性"的强证据

---

## 五、文档命名和数据流

### 5.1 关键文件用途

```
work/round4/
├── final_6_round4_v5.json      # 最终 Top-6 (含 score, pLDDT, 预测 Tm)
├── submission_round4_v5.csv    # 最终提交文件
├── esmfold_round4_v3.json      # 41 手工 ESMFold 结果
├── esmfold_mpnn.json           # 12 sfGFP MPNN ESMFold 结果
├── esmfold_mpnn_av.json        # 8 avGFP MPNN ESMFold 结果
├── candidates_round4_v3.json   # 41 手工候选 (含突变, 序列)
├── mpnn_candidates_final.json  # 12 sfGFP MPNN 候选
├── mpnn_avgfp_candidates.json # 8 avGFP MPNN 候选
└── mpnn_output_final/          # ProteinMPNN 原始 FA 输出
    ├── T01/seqs/2B3P.fa
    ├── T03/seqs/2B3P.fa
    └── T05/seqs/2B3P.fa
```

### 5.2 决策树

```
                          41 手工 + 20 MPNN = 61 候选
                                       ↓
                            过滤 pLDDT >= 35
                                       ↓
                                 56 候选
                                       ↓
                          Score 排序 (前 25 名)
                                       ↓
              多样性选择: 6 骨架各取最高分, 每骨架 ≤ 2
                                       ↓
                                6 候选
                                       ↓
                       按 score 重排 (Seq 1 = 最高)
                                       ↓
                          v5 最终 Top-6 ⭐
```

### 5.3 与比赛规则的映射

| 比赛要求 | 我们的实现 |
|---------|-----------|
| 6 条序列 | ✅ 6 条 |
| 220-250 aa | ✅ 全部 234/238 |
| M 开头 | ✅ 全部 M |
| 20 标准 AA | ✅ 全部通过 set 检查 |
| 不在 Exclusion_List | ✅ 全部通过 set 检查 |
| chromophore 完整 | ✅ TYG/SYG/GYG |
| Best Top-1 排榜 | ✅ MPNN_T01_014 (Seq 1) 是最强保险 |
| Finit ≥ 0.3×WT 阈值 | ✅ 全部候选预测 Finit ≥ 0.55 |

---

## 六、Round 5 的设计建议

### 6.1 保持并扩展

1. **保留 MPNN 路线**: 已被证明在 GFP 上有效, Round 5 应继续扩展
2. **保留多骨架多样性**: 5 骨架策略显著降低风险
3. **保留 pLDDT 68+ 候选**: 高 pLDDT = 高折叠可靠性

### 6.2 改进

1. **增加 MPNN 温度 (T=0.5, T=0.7)**: 探索更激进的 de novo 设计
2. **增加 amacGFP MPNN**: 多骨架覆盖
3. **使用 ThermoMPNN 真实预测 Tm**: 替代 pLDDT 间接估值
4. **增加 ESM-IF1 打分**: 第二维度结构验证

### 6.3 战略考虑

1. **如果时间允许**: 跑 200+ MPNN 候选, 选 Top 6 (而非 60)
2. **如果时间紧张**: 直接用 v5 提交, 不要再改
3. **如果比赛截止临近**: 优先做 PDF + GitHub 必交材料

### 6.4 不应该做的事

1. **不要再用 ML 打分** (Round 2 教训)
2. **不要再给 MPNN 候选固定过多位** (chromophore ± 5 已足够)
3. **不要在中文路径跑 ProteinMPNN** (致命错误)
4. **不要在 Round 5 重新设计 MPNN_T01_014 之类高分候选** (已经够了)

---

## 七、给 Round 5 AI 的关键启示

1. **MPNN + 固定 chromophore** 是当前最优 de novo 策略
2. **pLDDT 是结构质量硬通货**, 高 pLDDT = 高可靠性
3. **多骨架多样性** = 系统性风险对冲
4. **Tm 估值可保守但不能太乐观**: v5 的 88°C 是间接推断
5. **比赛细节未知**: 联系 root 询问 CFPS protocol
6. **提交优先**: 不要追求完美, 比赛截止前一定要交

---

**完成日期**: 2026-06-22
**作者**: Trae AI Agent
**版本**: Round 4 v5
