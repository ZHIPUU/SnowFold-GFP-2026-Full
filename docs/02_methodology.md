# 02 · 方法论(Methodology)

## 总体 Pipeline

```
┌─────────────────────────────────────────────────────────┐
│  Phase 0: 数据准备                                       │
│  - 141K GFP variant CFPS brightness 数据 (GFP_data.xlsx)│
│  - 5 个 WT 序列                                         │
│  - 13.5 万条排除列表                                     │
└─────────────────┬───────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────┐
│  Phase 1: 加性 Ridge 建模(per-parent)                   │
│  - 突变 → 单点 Δ brightness                             │
│  - R² per parent: 0.62-0.82                             │
│  - 组合搜索生成 97K 候选 → top-144                       │
└─────────────────┬───────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────┐
│  Phase 2-3: ESM2-150M PLL 评分 + Top-6                  │
│  - PLL log-likelihood 作为 fitness proxy                │
│  - 选 6 条提交 Round 1                                  │
└─────────────────┬───────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────┐
│  Phase 5-6: 设计思路 PDF + 综合分预测                   │
│  - 按比赛规则预测综合分                                  │
│  - Round 1 Seq 6 = 9.50(乐观)                          │
└─────────────────┬───────────────────────────────────────┘
                  │
                  ▼ (进入 Round 2)
┌─────────────────────────────────────────────────────────┐
│  Step 1: ESM2-650M 全量嵌入(141K × 1280)               │
│  - 下载 2.6 GB 模型                                     │
│  - 141K 序列 batch_size=128                             │
│  - 耗时 83 min,GPU mem 2.87 GB                          │
└─────────────────┬───────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────┐
│  Step 2: XGBoost GPU 训练 + 加性 Ridge baseline        │
│  - 1280-d ESM + type one-hot (4-d) → XGBoost GPU       │
│  - 9:1 stratified split                                 │
│  - val R²=0.5165,RMSE=0.661                            │
│  - 加性 Ridge baseline R²=0.222(各母体 R² 均为负)       │
└─────────────────┬───────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────┐
│  Step 3a: 现有数据天花板分析                              │
│  - 全库 XGBoost 预测                                    │
│  - max Finit/Finit_WT = 2.33× (avGFP)                   │
│  - 实际测量 max = 2.53×                                 │
│  - 结论: 数据天花板 <3×,必须叠加论文突变                 │
└─────────────────┬───────────────────────────────────────┘
                  │
                  ▼ (数据驱动天花板太低)
┌─────────────────────────────────────────────────────────┐
│  Step 3b (已废弃): Top-10K Ridge 候选重打分               │
│  - ESM2-650M 重嵌入 + XGBoost 打分                      │
│  - 结果: 所有 top 候选 finit_rel ≈ 0.000×                │
│  - 教训: 老 Ridge 搜索空间窄,被数据锁死                │
└─────────────────┬───────────────────────────────────────┘
                  │
                  ▼ (转向论文知识驱动)
┌─────────────────────────────────────────────────────────┐
│  Step 3c: 论文突变驱动手工设计 11 条候选                  │
│  - sfGFP 11 突变(亮度)                                 │
│  - TGP 7 稳定突变(Tm)                                  │
│  - TGP 5 表面 E 突变(抗聚集)                           │
│  - 4 个母体交叉组合                                     │
└─────────────────┬───────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────┐
│  Step 4: 轻量级验证(无 ESMFold)                         │
│  - chromophore 三联体检查(TYG/SYG/GYG)                 │
│  - 排除列表粗查(前 1000 条)                            │
│  - 论文突变计数评分                                     │
│  - 选 Top-6                                            │
└─────────────────┬───────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────┐
│  Step 5: 生成 Round 2 提交 CSV                          │
│  - 4 scaffold 覆盖(avGFP, amacGFP, sfGFP, ppluGFP,    │
│    cgreGFP)                                            │
│  - 全合规 + chromophore 完整                           │
└─────────────────┬───────────────────────────────────────┘
                  │
                  ▼ (Round 2 后续改进)
┌─────────────────────────────────────────────────────────┐
│  Step C: Epistasis 显式建模                              │
│  - 从 141K 提取 4595 个单点效应 + 48010 个双突变加性残差│
│  - 拼接特征: ESM + type + n_mut + ΣΔ + Σε + ...        │
│  - XGBoost GPU 重训 → val R² = 0.9162!                  │
│  - 接近往届顶尖 R²=0.9434                              │
└─────────────────┬───────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────┐
│  Step C2 (回退): 模型对 OOD 候选失效                     │
│  - 模型预测所有 sfGFP-style 候选 brightness <3.0       │
│  - 实际应 ≥4.5(基于论文先验)                           │
│  - 原因: 模型学到 "多突变 → 低亮度" 偏见                │
│  - 提交回退到 paper-knowledge Top-6                     │
└─────────────────────────────────────────────────────────┘
```

---

## 关键技术细节

### ESM2-650M 嵌入(Step 1)
- 模型:facebook/esm2_t33_650M_UR50D (33 层, 1280-d, 20 heads)
- 加载绕过 `torch.hub`(避免 HEAD 请求挂死)
- 直接 `torch.load(model_path)` + 手动构造 `ESM2(num_layers=33, embed_dim=1280, attention_heads=20, alphabet=Alphabet.from_architecture('ESM-1b'))`
- `strict=False` 容忍 ~570 missing keys(无影响)
- batch_size=128, ~28-30 seq/s, 83 min 总耗时
- GPU mem 2.87 GB(模型 2.6 GB + 嵌入 cache)
- 嵌入方式:`representations[33][:, 1:-1].mean(dim=1)`(去除 BOS/EOS,取残基平均)

### XGBoost GPU(Step 2)
- 参数:`tree_method='hist', device='cuda', max_depth=8, lr=0.05, subsample=0.85, colsample_bytree=0.7, min_child_weight=5, reg_alpha=0.1, reg_lambda=1.0`
- early_stop=50, num_boost_round=2000
- best_iter=1775 (84s 训练)
- 特征:1280 (ESM) + 4 (type one-hot) = 1284

### 加性 Ridge baseline(Step 2)
- per-parent 训练(4 个独立模型)
- 特征:0/1 突变存在性(sparse)
- Ridge alpha=1.0
- **关键发现**:per-parent R² 全为负(模型不如预测母体均值),但全局 R²=0.222(因为类型间差异大)

### Epistasis 特征工程(Step C)
```python
# 1. 单点效应估计
single_effects[(parent, from_aa, pos, to)] = mean(brightness - WT_baseline)
# 2. 双突变 epistasis 估计
pair_epistasis[(parent, frozenset([mut1, mut2]))] = mean(
    observed - WT_baseline - sum(single_effects)
)
# 3. 对每个变体,聚合特征:
# - single_sum: 所有单点效应之和
# - single_mean, single_max, single_min
# - pair_eps_sum: 所有对 epistasis 之和(查表)
# - pair_eps_count: 找到的 epistasis 对数
# - unseen_pair_count: 找不到 epistasis 的对数(模型用 ESM 处理)
```

### 论文突变(Step 3c 设计原则)
- **sfGFP 11 突变**(亮度):F64L/S65T/F99S/M153T/V163A + S30R/Y39N/N105T/Y145F/I171V/A206V
- **TGP 7 稳定突变**(Tm):K30I/A53S/T59P/V60A/T82A/K190E/K208R
- **TGP 5 表面 E 突变**(抗聚集):K45E/K73E/K117E/R149E/N158E
- 候选设计逻辑:不同母体 × 不同突变组合 = 跨 scaffold 多样性

---

## 与 Round 1 的关键差异

| 方面 | Round 1 | Round 2 |
|---|---|---|
| 嵌入模型 | ESM2-150M (480-d) | ESM2-650M (1280-d) |
| 训练架构 | 加性 Ridge | XGBoost GPU + epistasis 特征 |
| 候选搜索 | 数据驱动 (Ridge 组合) | **论文知识驱动** |
| scaffold 覆盖 | 5 个全用 | 5 个全用(更平衡) |
| 验证手段 | ESM PLL + 排除列表粗查 | chromophore + 排除列表 + 模型打分 |
| Tm 预测 | 粗略类别估计 | **无**(缺失能力) |
| pLDDT 验证 | 未做 | **未做**(ESMFold 不可用) |

---

## 失败的尝试(避免重蹈覆辙)

1. **Step 3b (Top-10K Ridge 重打分)**:彻底失败。所有 top 候选 finit_rel ≈ 0。原因:老 Ridge 搜索空间窄,只包含破坏折叠的突变堆砌。
2. **Step C2 (Epistasis 模型重打分)**:发现 OOD 陷阱。模型对训练分布内 R²=0.92 优秀,但对论文突变组合(分布外)系统性低估。
3. **ESMFold 加载**:fair-esm 2.0 不含, HuggingFace SSL 卡住,放弃。

---

*详细难点见 `docs/03_challenges.md`,文件清单见 `docs/08_appendix.md`*