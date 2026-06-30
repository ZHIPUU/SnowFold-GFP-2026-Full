# Round 4 下一步方向 — Next Steps

> **目标**:基于 Round 3 的成果和待解疑点,给出 Round 4 的优先级建议。先做高 ROI 的,再做锦上添花的。

---

## 🎯 立刻做的事(P0 - Round 4 第一天)

### S1. 重跑 ESMFold on GPU (10 分钟)

**问题**:Round 3 的 pLDDT 是 CPU 推理,绝对值偏低(40-50)。

**做法**:
```python
# 1. 确认 CUDA 可用
import torch
print(torch.cuda.is_available())  # 必须 True

# 2. 在有 CUDA 的 terminal 跑
python "D:\生信\2026Protein Design\work\round3\esmfold_validate_v2.py"
# (用 FP32 + GPU,准确 pLDDT)
```

**预期输出**:
- GPU 推理下 pLDDT 应该高 20-30 分
- 例如 sfGFP+I152S 可能从 48.3 → 70+
- 这会**改变最终候选排序**

**ROI**:高 — 10 分钟获得真实 pLDDT,可能完全改变 Round 4 候选集。

---

### S2. 联系 root 询问比赛细节 (5 分钟)

**问题**:Finit/Ffinal 的具体测量条件、Ffinal 加热时长、排名规则等均未知。

**做法**:直接问 root session。

**应问问题**:
1. Finit 在多低表达量下测量?
2. Ffinal 加热时长?温度?
3. 排名是 6 条平均还是 Best Top-1?
4. 排除列表是全序列匹配还是子串匹配?
5. Tm 是否直接测量,还是通过 Ffinal/Finit 推算?
6. 比赛截止日期?
7. Round 3 提交是否被接受?还是需要 Round 4?

**ROI**:极高 — 5 分钟获取的关键信息可能改变整个策略。

---

### S3. Tm 预测模型 (1-2 天)

**问题**:综合分公式中 Ffinal/Finit 没有直接预测能力。

**目标**:训练 ESM2-650M 嵌入 → Ridge 的 Tm 预测器。

**数据收集**:
- 文献中已知 GFP 的 Tm 值(sfGFP ~78°C, TGP ~85°C, StayGold ~92°C, mBaoJin ~92°C, mEGFP ~70°C, mCherry ~50°C 等)
- 至少 50-100 个 GFP 变体的 Tm 数据

**实施步骤**:
1. 收集 GFP Tm 数据(从 FPbase 数据库、Pédelacq 论文、Close 2015 论文等)
2. 用 ESM2-650M 嵌入每条序列(可复用 Round 2 的 esm650m_embeddings.npy)
3. Ridge 回归:Tm ~ ESM embedding
4. 验证集 R² 应该 > 0.7

**代码**:
```python
import numpy as np
from sklearn.linear_model import Ridge
from sklearn.model_selection import cross_val_score

# 加载嵌入
emb = np.load("work/round2/esm650m_embeddings.npy", mmap_mode="r")
# 加载 Tm 数据
tm_data = pd.read_csv("work/round3/tm_dataset.csv")  # 需要创建

X = emb[tm_data["emb_id"]]
y = tm_data["tm_celsius"]
model = Ridge(alpha=1.0)
scores = cross_val_score(model, X, y, cv=5, scoring="r2")
print(f"Tm R²: {scores.mean():.3f} ± {scores.std():.3f}")
```

**ROI**:高 — 综合分的另一半 Ffinal/Finit 没这个就不准。

---

## 🔧 短期改进(P1 - Round 4 第一周)

### S4. GeoEvoBuilder 生成 de novo 候选 (2-4 小时)

**问题**:Round 3 候选全部基于"已知突变组合",缺少 de novo 设计。

**工具**:GeoEvoBuilder (PNAS 2025, 北大来鲁华团队)
- 代码: github.com/PKUliujl/GeoEvoBuilder
- 已在 GFP 上验证有效(2.3× 荧光增强 + Tm 提升)
- 零样本设计,无需训练数据

**做法**:
1. 克隆 GeoEvoBuilder 仓库
2. 准备 PDB 结构输入(avGFP 2WUR, sfGFP 2B3P, cgreGFP 2HPW)
3. 固定 chromophore 三联体(65-67 残基)
4. 运行 GeoEvoBuilder 生成候选
5. ESMFold 验证
6. 排除列表检查

**预期收益**:
- 可能生成 Round 3 没有的全新候选
- 2.3× 增强的潜力

**风险**:
- 安装可能复杂(需要 PyRosetta 等依赖)
- 生成候选可能质量不一,需要筛选

**ROI**:中-高 — 2.3× 增强是显著优势。

---

### S5. 补回 mBaoJin 候选 (1-2 小时)

**问题**:Round 3 因 pLDDT 38.8 删除了 mBaoJin。

**做法**:
1. 尝试多个 mBaoJin 突变位点(E142D, V193I, L194M, D173N 等)
2. 对每个突变版本:
   - 全量排除列表检查
   - ESMFold pLDDT 验证
3. 选 pLDDT 最高的 mBaoJin 版本

**预期**:
- Tm ~92°C 是无与伦比的优势
- 即使 pLDDT 偏低,实际折叠仍可能正确

**替代工具**:
- 用 ColabFold API(可能对 StayGold 家族更准)
- 引用 PDB 8QBJ 的实测结构作为对照

**ROI**:中-高 — Tm 优势显著。

---

### S6. ProteinMPNN backbone 设计 (2-3 小时)

**问题**:Round 3 候选都是"已知突变组合",没尝试 de novo backbone 设计。

**工具**:ProteinMPNN (JACS 2024)
- 对自然蛋白 backbone 重新设计序列
- 可固定功能位点(65-67)保持活性
- 已有多个 GFP 应用案例

**做法**:
1. 准备 backbone PDB(avGFP 2WUR, sfGFP 2B3P)
2. 固定 chromophore 三联体
3. 运行 ProteinMPNN 生成 100+ 候选
4. ESMFold pLDDT 筛选(pLDDT > 80)
5. 排除列表检查

**预期收益**:
- 可能生成 Round 3 没有的全新结构
- ProteinMPNN 已验证可提高稳定性和表达

**ROI**:中 — 工具已成熟,但需要筛选。

---

### S7. 结构对齐 TGP → sfGFP (2-3 小时)

**问题**:TGP 突变基于 mAG 编号,与 sfGFP 编号不直接对应。

**做法**:
1. 用 PyMOL 或 Biopython 做结构对齐
2. mAG 3ADF vs sfGFP 2B3P vs avGFP 2WUR
3. 找出 TGP 7 稳定突变在 sfGFP 的等效位置
4. 在 sfGFP 基础上叠加 TGP 等效突变

**预期收益**:
- Tm 可达 ~85°C(若成功)
- 显著提升 Ffinal/Finit

**风险**:
- 跨家族结构对齐有不确定性
- 可能破坏折叠

**ROI**:中 — Tm 优势显著,但工作量大。

---

## 🌱 探索性(P2 - Round 4 第二周以后)

### S8. DCI_asym GNN — 动力学驱动 epistasis 预测 (3-5 天)

**工具**:Huynh et al. PNAS 2025
- 代码: github.com/SBOZKAN/GNN
- 用 DCI_asym(不对称动态耦合指数)替代查表式 epistasis
- 无需训练集 epistasis 标注

**优势**:
- 根本性解决 OOD 失效问题
- 可作为新 baseline

**风险**:
- 需要分子动力学预计算权重
- 实现复杂

**ROI**:高(长远),但投入大

---

### S9. 多任务学习(brightness + Tm + 折叠概率) (2-3 天)

**思路**:
- 共享 ESM2 embedding
- 3 个 XGBoost/MLP heads
- 数据:已知 GFP 的 (seq, brightness, Tm, fold_prob)

**收益**:
- Tm 预测能力(关键缺失)
- 折叠概率预测(可作为 pLDDT 替代)

**风险**:
- 需要收集足够的 Tm 数据
- 多任务训练调参复杂

**ROI**:高(长远)

---

### S10. ESM3 API(若可访问) (1 小时)

**工具**:EvolutionaryScale Forge API
- ESM3 多模态生成
- 可生成 de novo GFP

**限制**:
- 需要申请 API 访问
- 可能收费

**ROI**:中(若可访问)

---

### S11. 进化搜索 + 论文突变组合 (3-5 天)

**思路**:
- 用遗传算法或贝叶斯优化
- 起点:Round 3 的 6 条候选
- 变异:每代随机选 1-2 突变替换为 ESM 推荐的
- 适应度:pLDDT + chromophore + 排除列表

**风险**:
- 仍是 OOD 问题
- 但起点在 OOD 区域,变异可能扩散到更广空间

**ROI**:中

---

## 📊 决策矩阵(Round 4 选哪个?)

| 选项 | 估计耗时 | 估计 R² 提升 | 估计综合分提升 | 风险 | ROI |
|------|---------|------------|--------------|------|-----|
| S1 重跑 ESMFold on GPU | 10 min | 0 | 0(确认基础) | 极低 | ⭐⭐⭐⭐⭐ |
| S2 联系 root 询问细节 | 5 min | 0 | ?(可能改变策略) | 极低 | ⭐⭐⭐⭐⭐ |
| S3 Tm 预测模型 | 1-2 天 | 0 | 大(综合分更准) | 中 | ⭐⭐⭐⭐ |
| S4 GeoEvoBuilder | 2-4 h | 0 | 可能 +1-3 综合分 | 中 | ⭐⭐⭐⭐ |
| S5 补回 mBaoJin | 1-2 h | 0 | 可能 +0.5-1.5 综合分 | 中 | ⭐⭐⭐ |
| S6 ProteinMPNN | 2-3 h | 0 | 可能 +0.5-2 综合分 | 中 | ⭐⭐⭐ |
| S7 结构对齐 TGP → sfGFP | 2-3 h | 0 | 可能 +0.5-1 综合分 | 中 | ⭐⭐⭐ |
| S8 DCI_asym GNN | 3-5 天 | +0.02-0.05 | 可能 +0.5 综合分 | 高 | ⭐⭐ |
| S9 多任务学习 | 2-3 天 | +0.01 | 可能 +1 综合分 | 中 | ⭐⭐⭐ |
| S10 ESM3 API | 1 h | 0 | 可能 +1 综合分 | 低(若可访问) | ⭐⭐⭐ |
| S11 进化搜索 | 3-5 天 | 0 | 可能 +0.5 综合分 | 中 | ⭐⭐ |

---

## 🎯 推荐优先级

### 第一天(必须做)
1. **S1 重跑 ESMFold on GPU**(10 min)
2. **S2 联系 root 询问比赛细节**(5 min)

### 第一周(短期)
3. **S3 Tm 预测模型**(1-2 天,与其他并行)
4. **S5 补回 mBaoJin**(1-2 h)
5. **S7 结构对齐 TGP → sfGFP**(2-3 h)
6. **S6 ProteinMPNN**(2-3 h)

### 第二周及以后
7. **S4 GeoEvoBuilder**(2-4 h)
8. **S9 多任务学习**(2-3 天)
9. **S8 DCI_asym GNN**(3-5 天)

---

## 🔄 如果时间紧张(Round 4 截止前 1 天)

**最小可行方案**:
1. ✅ S1 重跑 ESMFold on GPU(10 min)
2. ✅ S2 联系 root 询问细节(5 min)
3. ✅ S5 补回 mBaoJin(1-2 h)
4. 根据 S1 结果决定:替换哪些候选?

**最大努力方案**:
1. ✅ S1, S2
2. ✅ S5
3. ✅ S7 结构对齐
4. ✅ S3 Tm 预测(粗略版)
5. 重新打分 + 选 6

---

## 📈 长期(3+ 月)

如果项目要继续,以下方向值得投入:

1. **AlphaFold2 集成 + 自定义 pLDDT 阈值**
2. **多任务学习(brightness + Tm + 折叠)**
3. **DCI_asym GNN 根本性解决 OOD**
4. **进化搜索(NSGA-II / CMA-ES)**
5. **多模态融合(MSA + 结构 + 序列)**
6. **强化学习做组合搜索**
7. **湿实验迭代(若可获得)**

但这些都是科研方向,**比赛不一定有时间投入**。

---

## 🎯 Round 4 候选集(预期)

如果 Round 4 完美执行,我们期望的最终候选集可能是:

| Seq | 候选 | 母体 | 突变数 | 预期 Tm | 预期亮度 |
|-----|------|------|--------|---------|---------|
| 1 | sfGFP + I152S | sfGFP | 1 | ~78°C | 中 |
| 2 | mBaoJin + E142D | mBaoJin | 1 | ~92°C ⭐ | 高 |
| 3 | avGFP + sfGFP 4 + TGP 4 等效 | avGFP | 8-10 | ~85°C ⭐ | 高 |
| 4 | GeoEvoBuilder de novo #1 | de novo | ? | 估 | 估 |
| 5 | ProteinMPNN design #1 | 设计 | ? | 估 | 估 |
| 6 | amacGFP + sfGFP 5 | amacGFP | 5 | ~75°C | 中 |

**预期综合分**:基于 Round 3 估计 4-7,Round 4 可达 6-9(加入 mBaoJin 的 92°C Tm 优势 + de novo 候选的潜力)。

---

*详见 [round3_05_open_questions.md](round3_05_open_questions.md) 了解待解疑点, [round3_07_handoff.md](round3_07_handoff.md) 了解接手指南。*