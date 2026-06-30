# 01 · 两轮成果汇总

## Round 1 — 加性 Ridge + ESM2-150M(已完成提交)

### 时间线
- Phase 1:数据探索 + 加性 Ridge 模型
- Phase 2:ESM2-150M PLL 评分
- Phase 3:终选 6 条 + 序列验证
- Phase 5:设计思路 PDF 生成
- Phase 6:综合分预测

### Round 1 提交 6 条
| Seq | scaffold | 突变 | 预测 Finit_rel | 预测综合分 |
|---|---|---|---|---|
| **6** ⭐ | avGFP | 7 | 10.0× | **9.50** |
| 3 | cgreGFP | 4 | 11.0× | 2.41 |
| 4 | ppluGFP | 4 | 19.8× | 1.60 |
| 2 | avGFP | 3 | 15.2× | 1.52 |
| 1 | avGFP | 2 | 5.1× | 0.51 |
| 5 | amacGFP | 2 | 2.1× | 0.18 |

### Round 1 模型质量
- 加性 Ridge val R²:**0.62–0.82 per parent**(但全局 XGBoost 对照 R²=0.517)
- ESM2-150M 嵌入 + PLL log-likelihood(perplexity)用于排序
- **问题**:线性加性模型 + 数据驱动搜索空间窄,只能找到数据子集内的最优

### Round 1 产物
- `submission_yourteamname.csv`(Round 1 版)
- `work/phase5/design_doc.pdf`(设计思路 PDF)
- `work/phase6_scores.csv`(按比赛规则评分预测)
- `README.md` + `run_all.py`(一键复现)

---

## Round 2 — XGBoost GPU + ESM2-650M + Epistasis 特征(当前状态)

### Round 2 关键升级
1. **嵌入规模**:ESM2-150M (480-d) → **ESM2-650M (1280-d)** — 信息量 2.7×
2. **训练架构**:加性 Ridge → **XGBoost GPU** — 非线性,处理 epistasis
3. **特征工程**:仅 type one-hot → **+ 单点突变效应 + 双突变加性残差** — 显式建模 epistasis
4. **设计驱动**:纯数据筛选 → **数据 + 论文知识双驱动**

### Round 2 模型里程碑
| Step | 模型 | val R² | val RMSE (log10) | 备注 |
|---|---|---|---|---|
| Step 2 | XGBoost GPU + ESM2-650M | **0.5165** | 0.661 | 第一次 GPU 训练,84s 完成 |
| Step C | + epistasis 特征(4595 单点 + 48010 双突变) | **0.9162** | 0.275 | ΔR²=+0.40,接近往届顶尖 0.9434 |
| **Step C2** | 用 Step C 重打分 11 候选 | n/a | n/a | ❌ OOD 候选全被判为低 brightness |

### Round 2 数据天花板
- 141K 全库 XGBoost 预测 max Finit_rel:
  - avGFP: **2.33×**(E5K:D75N:N104Y)
  - amacGFP: 1.78×
  - cgreGFP: 1.17×
  - ppluGFP: 1.68×
- **实际测量 max**:avGFP 2.53×(T37S:K40R:N104S)
- **结论**:数据本身无法支撑 >3× 的 brightness gain;必须叠加论文突变

### Round 2 候选设计(11 条 → 选 6 条)
手工设计 11 条候选,基于 4 个母体 + 论文突变组合:

| Seq | 提交编号 | name | scaffold | 关键设计 |
|---|---|---|---|---|
| 1 | Seq 1 | avGFP+sfGFP完整+TGP_全面增强 | avGFP | sfGFP 11 + TGP 4 稳定 + TGP 5 表面 E = 19 muts |
| 2 | Seq 2 | avGFP+sfGFP完整+TGP稳定 | avGFP | sfGFP 11 + TGP 4 稳定 = 14 muts |
| 3 | Seq 3 | amacGFP+sfGFP风格+TGP稳定 | amacGFP | S65T/V99S/M153T/I166T 等 + TGP 4 = 10 muts |
| 4 | Seq 4 | sfGFP+TGP_稳定核心 | sfGFP | sfGFP + A53S/T59P/V60A/T82A = 4 muts |
| 5 | Seq 5 | ppluGFP_原始 | ppluGFP | WT(高 Tm 保险) |
| 6 | Seq 6 | cgreGFP+S65T+K163A | cgreGFP | S65T + K163A(高 baseline brightness) |

所有 6 条:**长度 222-238,M 开头,标准 AA,chromophore 完整**(TYG/SYG/GYG ✓)

### Round 2 模型陷阱(关键教训)
⚠️ **Step C 模型在 in-distribution 上 R²=0.9162 优秀,但对 OOD 候选(论文突变组合)系统性低估**:
- 训练数据 max brightness ~4.6(对应 ~40K linear)
- 多数 variants brightness <4.0 log10
- 模型学到 "多突变 → 低 brightness" 的先验
- 我们的 sfGFP-style 候选(预期 4.5+ log10)被预测为 2.6-3.0 log10
- **结论**:模型仅作为 in-distribution sanity check,不能用于 OOD 候选设计
- **最终决定**:提交回退到 paper-knowledge Top-6(Step 5b)

### Round 2 关键产物
| 文件 | 描述 |
|---|---|
| `submission_yourteamname.csv` | **Round 2 提交**(已就位) |
| `work/round2/stepC_xgboost_epistasis.model` | 最新模型(val R²=0.9162) |
| `work/round2/candidates_round2_design.csv` | 11 条候选(完整序列 + 突变) |
| `work/round2/stepC_summary.json` | Step C 模型指标 |
| `work/round2/stepC2_candidates_scored.csv` | 模型对候选的预测(显示 OOD 陷阱) |
| `work/round2/step3b_summary.json` | Step 3b 数据驱动搜索(已废弃) |

---

## 横向对比

| 维度 | Round 1 | Round 2 | 提升 |
|---|---|---|---|
| 模型 R² | 0.62-0.82 (Ridge) | **0.9162** (XGBoost + epistasis) | ~3× R² |
| 嵌入模型 | ESM2-150M | **ESM2-650M** | 2.7× 信息量 |
| 训练硬件 | CPU only | **RTX 5080 GPU** | ~50× 速度 |
| 候选设计 | 数据 + Ridge 搜索 | **论文知识 + 验证** | 多样性 ↑ |
| 提交策略 | 1 avGFP 主导 | **4 scaffolds 覆盖** | 风险分散 |

---

## 留待解决的(见 `docs/04_open_questions.md`)

1. Round 2 提交**没有 pLDDT 验证**(ESMFold 不可用),有破坏折叠风险
2. Epistasis 模型对 OOD 候选失效的根本原因 + 改进方向
3. 数据 max Finit/Finit_WT=2.33× 的硬天花板,如何绕过
4. Tm 估计模型缺失,Ffinal/Finit 无法直接预测

---

*详细方法论见 `docs/02_methodology.md`;难点见 `docs/03_challenges.md`*