# Round 4 工作目录

## 比赛规则要点 (官方完整版, 2026-06-22 更新)

- **打榜**: 6 条提交序列中 **Best Top-1** 综合得分排名
- **综合分**: `Finitial/Finitial_WT × Ffinal/Finitial`
- **WT对照**: sfGFP (每批含对照)
- **极低亮度阈值**: Finitial < 0.3 × Finitial_WT (sfGFP) → 0 分淘汰
- **独立奖项**: 最佳亮度奖、最佳热稳定奖
- **序列约束**: 220-250 aa, M 开头, 仅 20 标准AA, 不在 Exclusion_List 中

## Round 4 战略调整 (Best Top-1 视角)

| 角色 | 数量 | 候选类型 | 目标 |
|------|------|---------|------|
| 🔥 热稳爆款 | 2 条 | mBaoJin 衍生 (Tm 92°C) | 综合分爆款 + 最佳热稳奖 |
| ⚖️ 综合平衡 | 2 条 | sfGFP + htFuncLib(sf:acid.3) | 综合分主力 |
| 🛡️ 保险条 | 1 条 | sfGFP + I152S (Round 3 验证) | 保底 |
| 🆕 探索 | 1 条 | sfGFP + TGP 等效突变 (结构对齐) | 高Tm尝试 |

## GPU 环境就绪状态

- ✅ torch 2.11.0+cu128
- ✅ RTX 5080 Laptop (16GB, sm_120 Blackwell)
- ✅ ESMFold 加载正常, GPU 推理 ~5s/序列 (CPU 130s)

## Round 3 → Round 4 关键修正

**Round 3 文档错误**: "CPU pLDDT 偏低 20-30 分, GPU 会高很多"

**GPU 实测**: CPU/GPU pLDDT 几乎一致 (差异 < 0.2)。pLDDT 40-50 是 ESMFold 对 GFP 家族的固有置信度水平, **不能用 50/70 作为硬阈值**, 只能用作相对排序参考。

## 脚本顺序

1. `01_seed_candidates.py` — 生成所有候选池 (基于论文知识)
2. `02_esmfold_screen.py` — GPU 批量 ESMFold 验证
3. `03_score_and_select.py` — 多指标综合打分, 选 Top-6
4. `04_finalize_submission.py` — 生成提交 CSV

## 关键文献依据

- **htFuncLib** (Weinstein 2023, Nat Commun): GFP active-site 设计, **Tm 达 96°C**, sf:acid.3 突变集已公开
- **mBaoJin** (Zhang 2024, Nat Methods): StayGold 单体, Tm 92°C, 8 突变 (S55T/H77R/E80G/Q140P/H141Q/C165Y/N171Y/T201A)
- **sfGFP** (Pédelacq 2006): 11 突变, Tm ~78°C, S30R +1.25 kcal/mol
- **TGP** (Close 2015): mAG 衍生, Tm ~85°C, 85°C 几乎完全保持荧光
- **Round 2 OOD 教训** (Ertelt 2025): ML 模型对论文突变组合系统性低估, **不用 ML 给 OOD 打分**
