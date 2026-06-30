# Round 6 总结报告

> **日期**: 2026-06-22
> **目标**: 重构打分体系，识别正确设计方向
> **状态**: ✅ 完成，识别出关键方法论错误并修正

---

## 一、本轮成果

### 1.1 核心认知突破

**Round 6 最大的收获不是新的提交序列，而是对设计方法论的根本性修正。**

| 修正前（错误） | 修正后（正确） | 来源 |
|---------------|---------------|------|
| pLDDT/pTM 作为"折扣因子" | 结构完整性是**硬性门槛** | 用户纠正 + 规则文件 |
| ML 模型（XGBoost）预测 OOD 亮度 | 结构门控优先，亮度后估 | 实验数据证实 OOD 失效 |
| mBaoJin 假设 1.8× 亮度 | pTM=0.37 → 结构崩塌 → 0 分 | ESMFold 验证 |
| MPNN 高突变数被惩罚 | 结构良好→不被低估，突变数无关 | 用户纠正 |
| 综合分可以很高（>1.0） | 受 pLDDT < 80 限制，合理范围 0.3-0.5 | 规则校准 |

### 1.2 修正后 Top-6 提交

[submission_round6_corrected.csv](file:///d:/生信/2026Protein%20Design/work/round6/submission_round6_corrected.csv)

| Seq | 名称 | 骨架 | pLDDT | 生色团 | pTM | 亮度 | Tm | **综合分** |
|:---:|------|:----:|:-----:|:-----:|:---:|:----:|:--:|:---------:|
| 1 | MPNN_T01_014 | sfGFP_MPNN | 68.3 | 67.6 | **0.765** | 0.44× | 81°C | **0.37** |
| 2 | G4_sfGFP_I152S_K166V_Q69L | sfGFP | 50.8 | 49.1 | 0.604 | 0.37× | 82°C | 0.34 |
| 3 | C1_sfGFP_I152S | sfGFP | 48.2 | 58.9 | 0.565 | 0.39× | 81°C | 0.33 |
| 4 | G1_sfGFP_I152S_Q69L_S72A | sfGFP | 48.6 | 54.0 | 0.574 | 0.37× | 81°C | 0.32 |
| 5 | MPNN_T01_024 | sfGFP_MPNN | 61.1 | **68.5** | 0.708 | 0.41× | 78°C | 0.30 |
| 6 | C2_sfGFP_I152S_Q69L | sfGFP | 47.9 | 50.7 | 0.567 | 0.34× | 81°C | 0.29 |

### 1.3 产出文件

| 文件 | 路径 |
|------|------|
| 结构差距分析 | [round7_gap_analysis.py](file:///d:/生信/2026Protein%20Design/work/round7_gap_analysis.py) |
| 修正版评分脚本 | [round6_final_corrected.py](file:///d:/生信/2026Protein%20Design/work/round6_final_corrected.py) |
| 结构分析脚本 | [round6_analyze_structures.py](file:///d:/生信/2026Protein%20Design/work/round6_analyze_structures.py) |
| 完整排名 | [round6_full_ranking.json](file:///d:/生信/2026Protein%20Design/work/round6/round6_full_ranking.json) |
| 最终 Top-6 | [final_6_round6_final.json](file:///d:/生信/2026Protein%20Design/work/round6/final_6_round6_final.json) |
| 提交 CSV | [submission_round6_corrected.csv](file:///d:/生信/2026Protein%20Design/work/round6/submission_round6_corrected.csv) |

---

## 二、与 Round 5 的差距

### 2.1 核心方法论改进

| 维度 | Round 5 | Round 6 |
|------|---------|---------|
| **评分逻辑** | 手工 Tm 估值 + ML 预测综合分 | 结构化门控 → 亮度估计 → 热稳估计 |
| **结构角色** | pLDDT 是众多因素之一 | pLDDT/pTM 是**硬性门槛** |
| **ML 模型使用** | XGBoost 预测绝对亮度（OOD 失效） | 仅用于已知序列区间，OOD 用结构推断 |
| **mBaoJin 处理** | 假设 1.8× 亮度, 92°C Tm | pTM<0.5 → 崩塌 → 0 分 |
| **MPNN 处理** | 57 突变得 0.14 分（被惩罚） | pTM=0.765 得 0.37 分（第一） |
| **阈值机制** | 无硬性结构门槛 | pTM<0.5=0分, pLDDT<80需改进 |

### 2.2 关键发现

1. **ML 模型（XGBoost）在 OOD 上完全不可靠**
   - 141K 训练数据 99% 为 0-3 个突变的 avGFP 变体
   - 对 sfGFP/MPNN 等多突变序列系统性预测亮度≈0
   - 结论：不能用 ML 模型预测训练分布之外的序列亮度

2. **104 条候选只有 6 条通过 pTM>0.5 的基本结构门控**
   - 通过数：pTM>0.75 (1), pTM>0.7 (6), pTM>0.5 (35), pTM<0.5 (69)
   - **0 条通过全部三项结构门槛**

3. **mBaoJin 全系结构崩塌**
   - 11 个 mBaoJin 候选 pTM 均 < 0.4（集中 0.36-0.38）
   - ESMFold 预测为无序状态 → 实验大概率测不到荧光

---

## 三、实现路径

### 3.1 Round 6 的执行流程

```
Step 1: 数据摸底
├── 读取 104 条去重候选 (来自 Round 4+5)
├── 读取 141K GFP_data.xlsx + esm650m_embeddings.npy
└── 确认训练数据以 avGFP 0-3 突变为主

Step 2: 用 XGBoost epistasis 模型评分
├── 发现 R²=0.916 但 OOD 完全失效
├── evolvepro_pred 对 sfGFP/MPNN 候选系统性低估
└── 确认 Round 2 就发现的 OOD 问题未解决

Step 3: 构建手工评分体系 v1
├── 文献亮度 × 结构折扣 × 热稳因子
├── mBaoJin 被高估 (1.35× → 1.29 分)
└── MPNN 被低估 (0.39× → 0.14 分)

Step 4: 用户纠正 → 重写 v2
├── 结构门控: pTM<0.5→0分
├── MPNN 不应被惩罚
├── 生色团质量决定亮度
└── 最终 Top-1: MPNN_T01_014 = 0.37 分

Step 5: 规则文件发布 → 精确定义门槛
├── pTM>0.75, pLDDT>80, 生色团>85 → 全部通过
├── 铜奖>0.30, 银奖>0.50, 金奖>0.80
└── 当前 0 条通过全部门槛
```

---

## 四、遇到的难点

### 4.1 OOD 亮度预测（未解决）

技术上最难的挑战。141K 数据集的训练分布决定了任何 ML 模型都无法准确预测 OOD 突变的亮度。这不是 tuning 能解决的——需要全新的训练数据或方法。

### 4.2 ESMFold 对非-avGFP 骨架的偏差

mBaoJin (StayGold 家族) 与 avGFP 序列同源性 <33%，ESMFold 训练数据中这类样本极少，导致系统性低 pLDDT。目前无简单解决方案。

### 4.3 esm v3 与 fair-esm v2 的包冲突

`esm` (v3.2.1) 和 `fair-esm` (v2.0.0) 同时安装，导致 `import esm` 导入错误的版本。ESMFold 无法通过 Python API 调用。

### 4.4 无候选通过全部三个结构门槛

当前 104 条候选的最佳生色团 pLDDT=68.5，距所需 85 还差 16.5。这在现有候选池中无法通过任何调参或重组解决——**必须重新设计**。

---

## 五、待解决的疑点

1. **ESMFold vs Boltz-1 的一致性**: 对设计蛋白，两者的 pLDDT/pTM 给出不同评分时，该信谁？（缺少实验验证参照）
2. **pLDDT 与荧光的相关性**: 对于 MPNN 设计的蛋白，pLDDT>80 是否真的意味着有荧光？目前纯属假设。
3. **生色团 pLDDT 的精准定义**: 残基 58-72 和 210-230 这个范围是否最优？ESMFold 的残基级 pLDDT 对这些位置的精度如何？
4. **ProteinMPNN 的 GFP β-barrel 适用性**: ProteinMPNN 训练数据中 β-barrel 的比例是否足够？之前低恢复率是否源于此？

---

## 六、下一轮方向

详见 Round 7 实施计划。核心思路：

> **"先确保结构，再追求功能"**

通过固定生色团残基 (Y66, G67, R96, E222) 的 ProteinMPNN 设计，预期将全局 pLDDT 从 68 提升至 80+，生色团 pLDDT 从 67 提升至 85+。

### 关键技术支撑
- **Fixed-residue ProteinMPNN**: 固定关键残基，只重设计其余部分（Sumida et al. JACS 2024）
- **EnhancedMPNN (ResiDPO)**: 用 pLDDT 作为 DPO 奖励信号，3× 设计成功率提升（ICLR 2026）
- **LigandMPNN**: 显式建模 chromophore 配体（Nat Methods 2025）
- **GeoEvoBuilder**: Zero-shot GFP 设计，晶体结构验证（PNAS 2025）

---

*Round 6 于 2026-06-22 完成。下一阶段重心从"评分配置"转向"序列生成"。*
