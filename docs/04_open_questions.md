# 04 · 待解疑点(Open Questions)

> 本文列出项目当前尚未解决、需要后续工作或决策的问题。按优先级排序。

---

## 🔴 P0 - 必须尽快解决(影响提交质量)

### Q1. Round 2 提交没有 pLDDT 验证

**状态**:Seq 1 (19 muts)、Seq 2 (14 muts) 等多条候选没有结构验证。chromophore 完整,但**整体折叠**未知。

**风险**:可能 Seq 1 折叠失败,根本无荧光。

**解决方向**:
- 离线下载 ESMFold 14GB 权重,手动绕过 SSL
- 或用本地 Rosetta / I-TASSER 预测结构
- 或退回到更保守的设计(只用 5-10 muts)

**建议**:至少对 Seq 1 做粗略结构检查,或换为更保守版本。

### Q2. Tm 预测能力缺失

**状态**:Ffinal/Finit 因子无法直接预测。我们用经验值(sfGFP ~78°C, StayGold ~92°C)估算,精度很差。

**影响**:综合分计算中,Ffinal/Finit 项基本是猜测。

**解决方向**:
- 训练专门的 Tm 预测模型(基于结构或序列)
- 用已知 Tm 的 GFP 库做迁移学习
- 引入 ΔΔG 预测(类似 FoldX 的 energy function)

**建议**:这是 Round 3 的核心改进方向。

---

## 🟡 P1 - 重要但非阻塞

### Q3. 数据天花板 (2.33×) 能否被绕过

**现状**:141K 数据 max Finit/Finit_WT=2.33×,根本达不到第一轮 9.50。

**根本问题**:数据基于"已部分 sfGFP 化的中间体",**真正的 sfGFP / TGP / StayGold 风格突变组合不在训练分布内**。

**可能的绕过**:
- **数据扩增**:用 ESM-IF1 等模型生成"虚拟组合",扩展训练集
- **迁移学习**:用文献报道的 sfGFP / TGP 亮度数据做 fine-tune
- **多任务学习**:同时预测 brightness + Tm + 折叠概率,共享 embedding

**现状**:用论文知识绕过了这个问题(手工设计候选),但模型没真正解决。

### Q4. Epistasis 模型 OOD 失效原因

**现象**:in-distribution val R²=0.92,但对论文候选预测全错。

**可能的根因分析**:
1. **训练分布偏度**:max brightness 4.60, mean ~3.5,模型学不到 ≥4.5 的区域
2. **特征冲突**:epistasis 特征是为训练集设计的(查表式),对 OOD 突变全是 0
3. **XGBoost 决策树偏好**:倾向于落在见过的特征空间

**TODO**:定量分析 — 在训练集中,log10 brightness 4.0+ 的样本占比多少?

### Q5. Step C2 提交为什么被回退

**决策记录**:Step C2 用 epistasis 模型重新打分 11 候选,选 Top-6 (Seq 11, 6, 1, 4, 5, 12)。但模型预测 finit_rel 全是 0.04-0.09×(应 ≥5×),**回退到 Step 5 的 paper-knowledge Top-6**。

**遗留问题**:如果用 Step C2 提交,实际综合分是多少?可能比 paper-knowledge 差,但**没有 ground truth**。

**建议**:在 Round 3 之前,**保留两种提交各一份**,作为对照实验。

---

## 🟢 P2 - 探索性,提升空间

### Q6. ESM2-3B / 15B 是否值得上

**信息**:往届顶尖选手 R²=0.9434, 我们 epistasis 模型 0.9162。差距 ~0.027。

**假说**:升级到 ESM2-3B (2560-d, 翻倍信息量)可能补上这 0.027 差距。

**风险**:
- ESM2-3B 模型 5.4 GB,需要重新嵌入 141K(预计 3-5 小时)
- 模型不一定能泛化到 OOD
- 我们的 6 条候选不在数据训练分布,3B 可能也只是稍微好一点

**建议**:有 4+ 小时空闲时再尝试。

### Q7. mBaoJin / StayGold 风格突变的精确清单

**现状**:TGP 有精确突变表(从 PDF Figure 2D alignment 提取),sfGFP 有 11 突变清单,**但 mBaoJin 和 StayGold 没有精确突变清单**。

**影响**:候选设计只能用 sfGFP + TGP 风格,不能叠加 mBaoJin / StayGold。

**解决**:
- 查 PDB 4TZA(TGP), 8A65(mBaoJin?), 8AB9(StayGold?)
- 做结构对齐,提取突变
- 或用 ESM-IF1 做 inverse folding

### Q8. 排除列表覆盖度

**现状**:只查了前 1000 条 / 13.5 万条(0.7%),可能漏检。

**解决**:全量检查(~10 min),确保 11 候选都不在排除列表。

**TODO**:写一个快速检查脚本。

### Q9. 跨母体突变兼容性

**疑问**:TGP 突变(mAG 编号)直接套用到 avGFP/cgreGFP 等,**位置是否真的对应?**

**风险**:不同 GFP 家族序列差异大,某些位置插入/缺失,直接套用可能破坏折叠。

**解决**:做结构对齐(avGFP vs mAG),确认 TGP 突变对应的 avGFP 位置。

---

## ⚫ 已决定但尚未验证

### D1. 提交策略:保守 vs 激进

**决策**:Round 2 选 4 母体覆盖 + 1 保守 ppluGFP。

**验证状态**:**未知**(等比赛结果)。

**回退方案**:Round 3 改用 6 条 avGFP 重设计(集中火力)。

### D2. Epistasis 模型只用于 in-distribution

**决策**:epistasis 模型仅作为 in-distribution sanity check,**不能用于 OOD 候选设计**。

**验证状态**:已通过 Step C2 反例验证(finit_rel 全 0.05-0.09×)。

---

## ⚫ 已放弃的方向

### X1. 纯数据驱动搜索

**放弃原因**:数据天花板 2.33×,无法达到 Round 1 9.50。

**证据**:Step 3b 重打分 10K 候选,所有 finit_rel ≈ 0.000×。

### X2. ESMFold 加载

**放弃原因**:fair-esm 2.0 不含 + HF SSL 问题 + 14GB 下载耗时。

**替代**:轻量级验证(chromophore + 排除列表 + 论文先验)。

---

## ❓ 信息缺失类问题

### M1. 比赛评分细节

- Finit/Finit_WT 在多低的 CFPS 表达量下测量?具体培养基/温度?
- Ffinal 加热时长?冷却条件?
- 排名是 6 条平均还是 Top-1?
- 排除列表的具体含义(全序列匹配还是子串)?

### M2. WT 真实表现

- 5 个 GFP WT 在比赛 CFPS 条件下的实际 brightness?
- 它们的 Tm 是多少?
- 排除列表里哪些是已知 high performers?

### M3. 上传渠道

- 提交 CSV 通过什么系统上传?
- 提交截止日期?
- Round 2 / Round 3 时间表?

### M4. 论文 PDF 完整突变表

- mBaoJin 论文 PDF 的精确突变清单?
- StayGold 论文 PDF 的精确突变清单?
- nature-Local fitness landscape PDF 里的高适应度 hotspot?

---

## 📋 待办清单(按 ROI 排序)

| 任务 | ROI | 估计时间 | 状态 |
|---|---|---|---|
| Q1: Round 2 pLDDT 验证 | 高(避免 0 分) | 2-3h | TODO |
| Q8: 全量排除列表检查 | 中 | 10 min | TODO |
| Q2: Tm 预测模型 | 高(影响综合分) | 1-2 天 | TODO |
| Q3: 数据扩增 | 中 | 4-8h | TODO |
| Q4: OOD 失效根因 | 中(理论) | 2h | TODO |
| Q6: ESM2-3B 升级 | 中-高 | 5h | TODO |
| Q7: mBaoJin / StayGold 突变 | 中 | 3-4h | TODO |
| Q9: 跨母体兼容性 | 中 | 1-2h | TODO |

---

*对应下一步方向见 `docs/05_next_steps.md`,对应论文知识见 `docs/06_paper_kb.md`*