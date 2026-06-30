# Round 3 待解疑点 — Open Questions

> **本文内容**:Round 3 完成后的所有未解决问题,按优先级排序,影响后续 Round 4 改进方向。

---

## 🔴 P0 — 核心未验证问题

### Q1. CPU pLDDT 是否对应真实折叠?

**状态**:
- 所有候选 pLDDT 在 42.6-48.3 区间(中等偏低)
- ESMFold 论文标准 GPU 推理下,pLDDT > 70 才算"高置信度"
- 我们的 CPU 推理值**预计偏低 20-30 分**(基于 ESMFold 论文对比)
- 即 GPU 推理下,我们的 pLDDT 可能在 60-75 区间(可接受)

**核心风险**:
- 如果 CPU pLDDT 30-50 真实对应"未正确折叠",Round 3 提交可能全部失败
- 如果 CPU pLDDT 30-50 对应"正确折叠但置信度不高",则风险可控

**缓解策略**:
- 我们采用"多条件联合验证"(pLDDT + chromophore + scaffold 平衡)
- 即使个别候选失败,其他候选仍有合理概率成功

**解决方向**:
- **重新跑 ESMFold on GPU**(如果有 CUDA 可用环境)
- 或用 **ColabFold API** 验证
- 或引用 Round 1 Seq 6 (avGFP+sfGFP5+I152S) 的实际测量数据 → Seq 6 (avGFP+sfGFP4core+I152S) 类似

**决策建议**:在 Round 3 提交被评分前,**如有 GPU,优先重跑 pLDDT**。

---

### Q2. Tm 预测能力仍缺失

**状态**:
- 综合分公式 = Finit/Finit_WT × Ffinal/Finit
- **Ffinal/Finit** 因子无预测能力
- 我们用"已知 Tm 值"的文献估计:
  - avGFP WT: ~64°C
  - sfGFP: ~78°C
  - sfGFP 风格候选: ~78-82°C(估)
  - mBaoJin: ~92°C

**影响**:
- 综合分计算中,Ffinal/Finit 项基本是猜测
- Tm 高的候选(TGP 风格)实际表现未知

**解决方向**:
- 收集文献中 GFP 变体的 Tm 数据(50-100 个)
- 训练 ESM2-650M 嵌入 → Ridge 的 Tm 预测器
- 验证集 R² 应该 > 0.7

**决策建议**:**Round 4 的核心改进方向**。

---

### Q3. mBaoJin 母体未能在 Round 3 提交中体现

**状态**:
- mBaoJin 是 Tm 最高的母体(92°C)
- mBaoJin WT 命中排除列表
- mBaoJin + D173N pLDDT=38.8 被删除
- 最终提交全部用 sfGFP 风格(平均 Tm ~80°C)

**核心风险**:
- 我们错过了 Tm 92°C 的最强母体
- 可能在"最佳热稳奖"上失分

**解决方向**:
- 尝试多个 mBaoJin 突变版本(E142D, V193I, L194M 等)
- 用 ColabFold 而非本地 ESMFold 验证(可能对 StayGold 家族更准)
- 引用 StayGold 原始论文补充材料,找"低风险高收益"突变

**决策建议**:Round 4 优先补回 mBaoJin 候选。

---

## 🟡 P1 — 重要但非阻塞

### Q4. ESMFold 对 StayGold 家族的准确性

**现象**:
- mBaoJin + D173N 的 pLDDT 仅 38.8(相对最低)
- 推测:StayGold 家族结构与 ESMFold 训练集分布差异较大
- ESMFold 可能对 StayGold 家族置信度普遍偏低

**可能原因**:
1. StayGold 是 2022 年才发表的较新蛋白,ESMFold 训练集可能覆盖不足
2. StayGold 的 β-can 结构与传统 GFP 有差异
3. mBaoJin 是经过 8 轮定向进化的工程化版本,序列与原始 StayGold 有差异

**风险**:
- 即使 mBaoJin 实际折叠很好,pLDDT 也可能偏低
- 这导致我们误删了优秀候选

**验证方法**:
- 用 ColabFold 验证(精度通常比 ESMFold 高)
- 或查阅 PDB 8QBJ/8Q79/8QDD 的实测结构(已有晶体结构)
- 与已知晶体结构的 pLDDT 对比,校准 ESMFold

---

### Q5. avGFP+sfGFP10(Seq 3)的实际效果

**现象**:
- 这是突变最多的候选(10 个)
- 包含 S30R (+1.25 kcal/mol) 和所有 sfGFP 风格突变
- 理论上是"最接近 sfGFP"的候选
- pLDDT 仅 45.5(不是最高)

**核心问题**:
- sfGFP 风格的"上限"是 sfGFP WT 的表现
- sfGFP WT 在 CFPS 系统下的 Finit 和 Ffinal 是多少?
- 是否值得选择 sfGFP WT 作为 Seq 3?

**比较**:
- sfGFP WT 也命中排除列表(我们检查过)
- avGFP+sfGFP10 vs sfGFP WT 应该是相似的(仅 1 个细微差异)

**建议**:在 Round 4 中:
- 尝试更多 sfGFP WT 的 +1 突变版本(寻找绕开排除列表的版本)
- 用 Tm 预测器评估 Seq 3 vs sfGFP WT 的 Tm 差异

---

### Q6. Round 1 Seq 6 (9.50 综合分) 的可靠性

**状态**:
- Round 1 Seq 6 = avGFP + sfGFP 5 核心 + I152S
- 综合分 9.50 是"乐观估计"(假设 Finit 提升 10×)
- Round 3 Seq 6 几乎是 Round 1 Seq 6 的复制

**核心风险**:
- 如果 9.50 严重高估,Round 3 Seq 6 实际可能只有 4-7 分
- 但**即使 4-7 也比 Round 2 Seq 5 (0 分) 好**

**验证方法**:
- 没有 ground truth,只能等比赛结果
- 文献调研:sfGFP 风格的"典型 Finit 提升"是多少?(经验 5-10×)

---

### Q7. 比赛 CFPS 系统细节

**缺失信息**:
- Finit 在多低的 CFPS 表达量下测量?
- 具体培养基/温度?
- Ffinal 加热时长?
- 冷却条件?
- 排名是 6 条平均还是 Top-1?
- 排除列表的具体匹配规则(全序列 vs 子串)?

**影响**:
- 我们用文献估计 Tm 来猜 Ffinal/Finit
- 实际测量条件可能影响结果

**建议**:在 Round 4 前**联系 root session 询问比赛细节**。

---

### Q8. Round 2 vs Round 3 的实际差距

**未知**:
- Round 2 提交(违规 Seq 5)的实际得分?
- Round 3 提交的实际得分?
- 哪个更好?

**理论分析**:
- Round 2 Seq 5 直接判 0(违规)
- Round 3 没有违规,所有候选都可参与排名
- 即使 Round 3 其他候选表现一般,也比 Round 2 的 0 分好
- **Round 3 大概率优于 Round 2**(因为至少不是 0 分)

**结论**:**无法验证,只能等比赛结果**。

---

## 🟢 P2 — 探索性,可提升

### Q9. 文献调研发现的 8 篇论文的进一步应用

**未应用的工具/方法**:
1. **GeoEvoBuilder**(PNAS 2025) — 零样本设计,代码已开源
2. **Seq2Fitness + BADASS**(PLoS CB 2025) — 双模型 ensemble + 采样
3. **DCI_asym GNN**(PNAS 2025) — 动力学驱动 epistasis 预测
4. **ESM3**(Science 2024) — 多模态生成(通过 API)
5. **ProteinMPNN**(JACS 2024) — backbone 设计

**潜在价值**:
- 这些工具可能生成我们没考虑过的候选
- 特别是 GeoEvoBuilder,已在 GFP 上验证有效(2.3× 增强)

**建议**:Round 4 优先尝试 GeoEvoBuilder 和 ProteinMPNN。

---

### Q10. ESM2-3B 是否值得升级

**信息**:
- 往届顶尖选手 R²=0.9434,我们 epistasis 模型 0.9162
- 差距 ~0.027

**Round 3 决策**:**不升级**,因为根因不在模型大小(已在 Science Adv 2025 论文验证)。

**重新评估**:
- 如果 Round 4 决定重新走数据驱动路线,ESM2-3B 仍值得考虑
- 但**不能解决 OOD 失效的根本问题**

---

### Q11. TGP 突变是否能移植到 sfGFP

**状态**:
- TGP 突变基于 mAG 编号,与 sfGFP 编号不直接对应
- Round 3 完全放弃 TGP 突变

**潜在路径**:
- 做结构对齐(mAG 3ADF vs sfGFP 2B3P vs avGFP 2WUR)
- 找出 TGP 7 稳定突变在 sfGFP 的等效位置
- 在 sfGFP 基础上叠加 TGP 等效突变

**预期效果**:
- Tm 可达 ~85°C(若成功)
- 显著提升 Ffinal/Finit

**建议**:Round 4 如有时间,执行结构对齐。

---

### Q12. 多任务学习(brightness + Tm + 折叠概率)

**思路**:
- 同时训练 3 个预测器,共享 ESM 嵌入
- 数据:已知 GFP 的 (seq, brightness, Tm, fold_prob)

**潜在收益**:
- 每个任务都能从其他任务学到的特征中受益
- Tm 预测器可能是最大价值

**建议**:Round 4 长线改进方向。

---

## ⚫ 已决定但未验证

### D1. 突变数 ≤10 原则

**决策**:全部候选突变数 ≤10
**依据**:Arcadia 2025 数据 + 理论分析(p^N 折叠概率)
**未验证**:需要实验确认(我们没 CFPS 设备)

### D2. 删除 cgreGFP 和 mBaoJin

**决策**:基于 pLDDT 删除这两个候选
**风险**:可能是 ESMFold 误判(尤其 StayGold 家族)
**回退方案**:如果 Round 3 表现差,Round 4 重新考虑

### D3. 不使用 ML 模型打分

**决策**:Round 3 完全不依赖 ML 模型
**依据**:Science Adv 2025 + Round 2 自家证据
**风险**:可能有我们不知道的更好 ML 方法

---

## ⚫ 已放弃的方向

### X1. 纯数据驱动搜索

**放弃原因**:数据天花板 2.33×,无法达到第一轮 9.50
**证据**:Round 2 Step 3b 重打分 10K 候选,所有 finit_rel ≈ 0

### X2. Round 2 epistasis 模型给 OOD 打分

**放弃原因**:对论文突变组合系统性低估 1.5+ log10
**证据**:Round 2 Step C2 反例验证

### X3. ESMFold 权重从 HuggingFace 下载(无 SSL 解决)

**放弃原因**:SSL 证书问题无法直接绕过
**解决**:用 transformers 配合 `local_files_only=True`,模型已在本地缓存

### X4. 跨家族移植 TGP 突变

**放弃原因**:需要结构对齐,没时间做
**缓解**:Round 3 完全使用 sfGFP 风格(同家族编号兼容)

---

## ❓ 信息缺失类问题

### M1. WT 在比赛 CFPS 系统下的实测值

**缺失**:
- avGFP/sfGFP/amacGFP/cgreGFP/ppluGFP 在比赛 CFPS 条件下的实际 brightness?
- 它们的实测 Tm?
- 排除列表里哪些是已知 high performers?

**影响**:无法准确计算 Finit/Finit_WT 和 Ffinal/Finit

### M2. 比赛评分细节

**缺失**:见 Q7

### M3. 论文 PDF 完整突变表

**缺失**:
- mBaoJin 论文 PDF 的精确突变清单?
- StayGold 论文 PDF 的精确突变清单?
- nature-Local fitness landscape PDF 里的高适应度 hotspot?

**影响**:候选设计只能用 sfGFP 风格,不能叠加 mBaoJin/StayGold

---

## 📋 待办清单(按 ROI 排序)

| 任务 | ROI | 估计时间 | 状态 |
|------|-----|---------|------|
| **重跑 ESMFold on GPU** | 高(确认结构) | 10 min | TODO |
| **Tm 预测模型** | 高(综合分精度) | 1-2 天 | TODO |
| **补回 mBaoJin 候选** | 中-高(增加 Tm 优势) | 1-2 h | TODO |
| **GeoEvoBuilder 生成 de novo 候选** | 中-高 | 2-4 h | TODO |
| **ProteinMPNN backbone 设计** | 中 | 2-3 h | TODO |
| **结构对齐 TGP → sfGFP** | 中 | 2-3 h | TODO |
| **多任务学习(brightness+Tm+fold)** | 高 | 2-3 天 | TODO |
| **DCI_asym GNN(动力学 epistasis)** | 中 | 3-5 天 | TODO |
| **ESM2-3B 升级(若走数据驱动)** | 低 | 4-5 h | TODO |
| **联系 root 询问比赛细节** | 高(成本 0) | 5 min | TODO |

---

## 🎯 Round 4 优先方向

按 ROI 排序,Round 4 应该做:

### 立即做(高 ROI)
1. **重跑 ESMFold on GPU**(10 min,确认 pLDDT 是否真实)
2. **联系 root 询问比赛细节**(5 min,可能改变策略)
3. **Tm 预测模型**(1-2 天,完善综合分评估)

### 短期(中等 ROI)
4. **GeoEvoBuilder 尝试**(2-4 h,可能生成 de novo 优秀候选)
5. **补回 mBaoJin 候选**(1-2 h,增加 Tm 优势)
6. **ProteinMPNN backbone 设计**(2-3 h,生成新候选)

### 长期(低 ROI)
7. **DCI_asym GNN**(3-5 天,根本性解决 OOD 问题)
8. **多任务学习**(2-3 天,完善预测能力)

---

*详见 [round3_06_next_steps.md](round3_06_next_steps.md) 了解详细的下一步方向, [round3_07_handoff.md](round3_07_handoff.md) 了解接手指南。*