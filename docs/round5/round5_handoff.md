# Round 5 交接文档

> **日期**: 2026-06-22
> **状态**: Round 5 完成，未产生突破性进展
> **核心结论**: LigandMPNN/ProteinMPNN 工具链已搭建，但**未生成比 Round 4 更好的候选**。当前最佳候选仍是 Round 4 的 MPNN_T01_014 (pLDDT 68.3, pTM 0.765)。
> **给下一个 Agent 的警告**: 请先读"自我批评"部分，避免重复我们的错误。

---

## 一、Round 5 做了什么

### 1.1 文献调研（完成）

检索了 33+ 篇 2024-2026 前沿论文，整理了完整文献库：

- **文献库**: `work/round5/literature/round5_literature_kb.md`
- **需订阅下载清单**: `work/round5/literature/PAPERS_TO_DOWNLOAD.md`
- **已下载 4 篇 PDF**: `work/round5/文献/`（ESM3、GeoEvoBuilder、LigandMPNN、ESM3 SI）

关键发现：
- **LigandMPNN** (Nature Methods 2025) — 把 chromophore 当配体感知，比 ProteinMPNN 更适合 GFP
- **EVOLVEpro** (Science 2024) — ESM-2 + RF 顶层回归器，few-shot 定向进化
- **GeoEvoBuilder** (PNAS 2025) — 零样本同时优化稳定+活性，GFP 上已验证
- **htFuncLib** (Nat Commun 2023) — 16K 功能性 GFP 变体，Tm 最高 96°C
- **sfGFP 真实 Tm = 86.1°C**（BMC Res Notes 2023 CD 测定），不是之前以为的 78°C

### 1.2 LigandMPNN 部署（完成）

- 克隆了 LigandMPNN 仓库 + 下载了 4 个模型权重
- 修复了 `np.int` 兼容性问题（Python 3.14 + numpy 2.x）
- 在 avGFP (2WUR) 上成功运行：**190 条候选，chromophore 100% 保留**
- 多温度多种子扩展：T=0.1/0.15/0.2/0.25/0.4 × seed=37/137/237

### 1.3 EVOLVEpro 风格 fine-tune（完成）

- 用比赛 141K 数据训练了 XGBoost GPU 顶层回归器（ESM-2 650M 嵌入）
- val R² ≈ 0.5（与 Round 2 baseline 一致，不如 epistasis 模型的 0.916）
- 对 104 条候选做了独立亮度预测
- **结论**: EVOLVEpro 预测值全部偏低 (0.1-0.4 log10)，OOD 问题依然存在

### 1.4 ESMFold 批量验证（完成）

- 对 104 条候选全部做了 GPU ESMFold 评估
- 新增 LigandMPNN 候选最高 pLDDT = 64.8 (R5e_av_lmpnn_T015_s237_015, pTM 0.743)
- **但未超过 Round 4 的 MPNN_T01_014** (pLDDT 68.3, pTM 0.765)

### 1.5 GeoEvoBuilder（未完成）

- 克隆了仓库，但**模型权重在北大网盘**（访问码 `xx7W`），需要手动下载 `Se.pt`
- 代码需要 conda 环境 + dssp + libboost，当前环境不满足

---

## 二、对 Round 4 的改进

| 维度 | Round 4 | Round 5 | 改进 |
|------|---------|---------|------|
| 候选池规模 | 61 | 104 | +70% |
| 设计工具 | ProteinMPNN | + LigandMPNN (chromo-aware) | 新工具 |
| 亮度预测 | 无 | EVOLVEpro XGBoost (val R²~0.5) | 新维度 |
| 文献覆盖 | 5 篇参考论文 | 33+ 篇 2024-2026 前沿 | 大幅扩展 |
| 最高 pLDDT | 68.3 (MPNN_T01_014) | 68.3 (未超越) | **持平** |
| 最高 pTM | 0.765 | 0.765 (未超越) | **持平** |
| 提交合规 | ✅ 6/6 | ✅ 6/6 | 持平 |

**诚实评价**: Round 5 在工具链和文献上做了大量基础设施工作，但**没有生成比 Round 4 更好的候选序列**。LigandMPNN 虽然正确保留了 chromophore，但其最高 pLDDT (64.8) 仍低于 ProteinMPNN (68.3)。

---

## 三、当前最佳提交

### 混合策略 Top-6

| Seq | 名称 | 骨架 | 突变 | pLDDT | pTM | 角色 |
|-----|------|------|------|-------|-----|------|
| 1 | MPNN_T01_014 | sfGFP_MPNN | 57 | 68.3 | 0.765 | 冲刺 |
| 2 | R5e_av_lmpnn_T015_s237_015 | avGFP_LMPNN | 173 | 64.8 | 0.743 | 冲刺 |
| 3 | G4_sfGFP_I152S_K166V_Q69L | sfGFP | 3 | 50.8 | 0.604 | 保底 |
| 4 | X4_avGFP_sfGFP_I152S_Q69L | avGFP | 12 | 49.8 | 0.592 | 保底 |
| 5 | G1_sfGFP_I152S_Q69L_S72A | sfGFP | 3 | 48.6 | 0.574 | 保底 |
| 6 | Z3_amacGFP_sfGFP5_I152S | amacGFP | 5 | 45.5 | 0.534 | 保底 |

提交文件: `work/round5/submission_round5_hybrid.csv`

---

## 四、自我批评（给下一个 Agent 的警告）

### 问题 1: 缺乏新进展（循环回收）

**事实**: Seq 1 和 Seq 2 与 Round 4 完全相同。Seq 3-6 来自 Round 4 的手工设计候选池。Round 5 没有生成任何"新的更好序列"。

**根因**: LigandMPNN 生成了 190 条新候选，但 ESMFold 验证后最高 pLDDT 仅 64.8，不如 Round 4 的 ProteinMPNN (68.3)。新工具没有带来质量提升。

**教训**: 工具升级 ≠ 结果升级。LigandMPNN 在 avGFP 2WUR 上跑得好（chromo 100% 保留），但生成的序列质量不如 ProteinMPNN 在 sfGFP 2B3P 上的结果。

### 问题 2: 无效的组合诱变（拼凑的垃圾）

**事实**: Seq 3-5 引入了 I152S + Q69L + K166V + S72A 等突变组合，这些突变各自可能有益，但**组合后 pTM (0.57-0.60) 甚至比 Round 4 的 sfGFP+I152S (pTM 0.565) 没有显著提升**。

**根因**: 这些突变来自不同论文、不同骨架、不同设计目标，直接组合破坏了上位协同网络 (epistatic network)。

**教训**: 不要把不同来源的突变简单叠加。Round 2 的教训"突变堆砌"在 Round 5 重演了。

### 问题 3: 结构阈值依然失败

**事实**: 除 Seq 1 (pTM 0.765) 外，其余 5 条 pTM 均 < 0.65。pLDDT 均 < 70，侧链高度不确定。

**根因**: 手工设计的低突变候选 (3-12 mut) 虽然"安全"（大概率通过 30% 亮度阈值），但 ESMFold 对它们的结构信心不高。这说明 sfGFP/avGFP 家族在 ESMFold 中本身 pLDDT 就偏低（~45-50），不是因为设计不好。

**教训**: ESMFold 对 GFP 家族的 pLDDT 普遍偏低（WT sfGFP 本身也只有 ~50），不应该用 pLDDT 70 作为 GFP 的"正确折叠"门槛。**pTM > 0.5 可能是更合理的门槛**。

### 问题 4: 策略停滞

**事实**: Round 5 仍然在手工/半理性调整已知点突变（I152S, Q69L, S72A 等），没有做真正的 constrained MPNN de novo 设计。

**根因**:
- LigandMPNN 在 sfGFP (2B3P) 上失败（CRO chromophore 被 MPNN 重设计）
- LigandMPNN 在 avGFP (2WUR) 上成功但 pLDDT 不如 ProteinMPNN
- GeoEvoBuilder 缺权重文件未跑
- ESM3 API 未申请
- RFdiffusion3 未部署

**教训**: 工具部署耗时远超预期（np.int bug、中文路径 bug、模型下载、依赖冲突），导致实际设计时间不足。

---

## 五、实现路径

### 5.1 环境搭建

```
Python 3.14 + PyTorch 2.11.0+cu128 (RTX 5080, sm_120)
- ESMFold: transformers 4.48.1 (GPU FP32, 5s/序列)
- ProteinMPNN: 已部署 (Round 4)
- LigandMPNN: 已部署 (Round 5, 修复 np.int + ProDy 依赖)
- XGBoost GPU: 已部署 (Round 2, 复用)
- ESM-2 650M: 已下载 (hf-mirror.com)
```

### 5.2 候选生成 Pipeline

```
Round 5 完整 Pipeline:
1. PDB 下载 (2B3P, 2WUR, 1GFL) → RCSB
2. LigandMPNN 设计 (avGFP 2WUR, chromo 固定)
   → 190 条候选, 100% chromo 保留
3. ESMFold GPU 验证 (30 条 top, ~5s/条)
   → 最高 pLDDT 64.8, pTM 0.743
4. EVOLVEpro fine-tune (141K 数据 → XGBoost GPU)
   → val R² ~0.5, OOD 预测偏低
5. 多维度评分 (pLDDT + pTM + EVOLVEpro + LigandMPNN confidence)
6. Top-6 选择 (2 冲刺 + 4 保底)
```

### 5.3 关键代码文件

| 文件 | 用途 |
|------|------|
| `01_ligandmpnn_design.py` | LigandMPNN 初版设计 (CRO 问题) |
| `03_ligandmpnn_v2.py` | LigandMPNN 修正版 (chromo 固定) |
| `08_lmpnn_expanded.py` | LigandMPNN 多温度扩展 (190 条) |
| `09b_evolvepro_xgb.py` | EVOLVEpro XGBoost 训练 |
| `13_evolvepro_score.py` | EVOLVEpro 候选评分 |
| `17_corrected_selection.py` | 修正版 Top-6 选择 (pTM 优先) |

---

## 六、遇到的难点

### 6.1 LigandMPNN CRO 问题
- **问题**: 2B3P 的 chromophore 是 CRO (HETATM)，LigandMPNN 把它当可设计残基全部重写
- **解决**: 改用 2WUR (avGFP, 含原生 TYG 三联体)，固定 pos 30-42
- **教训**: LigandMPNN 对 HETATM 残基的处理与预期不同

### 6.2 numpy 2.x 兼容性
- **问题**: LigandMPNN 的 openfold 代码用 `np.int`，numpy 2.x 已移除
- **解决**: 批量替换 `np.int` → `int`
- **教训**: Python 3.14 + numpy 2.x 对老代码兼容性差

### 6.3 ESM-2 650M 模型下载
- **问题**: HuggingFace 默认源不通，`local_files_only=True` 失败
- **解决**: 设置 `HF_ENDPOINT=https://hf-mirror.com` 从镜像下载
- **教训**: 国内环境需要配置 HF 镜像

### 6.4 EVOLVEpro OOD 问题
- **问题**: XGBoost 训练 val R²~0.5，但对候选预测值全部偏低 (0.1-0.4 log10)
- **根因**: 候选远离训练分布（141K 数据的突变数 ≤5，我们的候选 3-173 突变）
- **教训**: 与 Round 2 epistasis 模型 OOD 失败一致，ML 在 OOD 上无法给出可靠亮度预测

---

## 七、待解决的疑点

### Q1: MPNN_T01_014 的 30% 阈值风险
- 57 个突变，pLDDT 68.3（高），但亮度可能 <30% WT
- **没有可靠方法预测**（EVOLVEpro 给出 0.12 log10，远低于训练集均值 3.3）
- 如果它折叠正确但 chromophore 不成熟 → 亮度低 → 0 分

### Q2: LigandMPNN vs ProteinMPNN 谁更好
- LigandMPNN pLDDT 最高 64.8 < ProteinMPNN 68.3
- 但 LigandMPNN 有 chromophore-aware 评分（ligand_confidence）
- **未做实验验证，无法确定哪个在 CFPS 中表现更好**

### Q3: ESMFold 对 GFP 的 pLDDT 偏低
- sfGFP WT 本身 pLDDT ~48，pTM ~0.55
- AlphaFold2/ColabFold 可能给出更高分数（有 MSA 信息）
- **未用 ColabFold/Boltz-1 交叉验证**

### Q4: mBaoJin 是否真的不可用
- pTM 0.36 看起来是"无序卷曲"
- 但 PDB 8QBJ 晶体结构已验证 mBaoJin 正确折叠
- **可能是 ESMFold 对 StayGold 家族的系统性偏差**
- Boltz-1 或 ColabFold 可能给出不同结论

### Q5: htFuncLib 16K 数据是否可用
- Addgene 有 plasmid，论文 SI 有序列表
- 16K 设计含实测 Tm（最高 96°C）和亮度
- **未下载和使用这批数据做 fine-tune**

### Q6: 比赛实际 CFPS 条件
- 加热时长未知（假设 30 min）
- CFPS 配方未知
- 排名是 Best Top-1 还是 6 条平均（规则已确认是 Best Top-1）
- **加热时长对 Ffinal/Finit 影响极大**

---

## 八、下一步方向（按 ROI 排序）

### P0: 必须做（可能产生突破）

1. **GeoEvoBuilder 复刻 1GFL 设计**
   - 下载权重 `Se.pt`（北大网盘，访问码 `xx7W`）
   - 论文已验证 GFP 设计有效（1GFL-15, 1GFL-19 实验确认荧光+热稳）
   - 这是**唯一在 GFP 上实验验证过的 zero-shot 设计工具**
   - **预期**: 可能生成全新的高质量候选

2. **ColabFold / Boltz-1 交叉验证**
   - 对 MPNN_T01_014 和 mBaoJin 候选做 ColabFold 预测
   - ColabFold 用 MSA 信息，可能给出比 ESMFold 更高的 pLDDT
   - 如果 mBaoJin 在 ColabFold 中 pTM > 0.5，说明 ESMFold 有偏见

3. **htFuncLib 16K 数据下载 + fine-tune**
   - 从论文 SI 或 Addgene 获取 16K 设计序列 + 实测 Tm
   - 用这批数据 fine-tune ESM-2（EVOLVEpro 风格）
   - **可能解决 OOD 问题**（htFuncLib 数据在 GFP 高突变空间）

### P1: 重要（稳健提升）

4. **ProteinMPNN 在更多 PDB 上跑**
   - 1GFL (GFP WT 原始晶体) — 可能比 2B3P (sfGFP) 给出不同结果
   - 8QBJ (mBaoJin) — 但 chromophore 格式问题
   - 7LG4 (amacGFP) — 跨骨架多样性

5. **ESM3 Forge API esmGFP 设计**
   - 申请 Forge API token
   - 复刻 Hayes 2024 的 chain-of-thought GFP 生成
   - **可能生成 58% 同源性的全新荧光蛋白**

6. **多目标 Pareto 优化 (MOG-DFM / NSGA-II)**
   - 同时优化 pTM + 突变数 + EVOLVEpro + Tm
   - 给 104 条候选做 Pareto 排序

### P2: 探索性

7. **RFdiffusion3 all-atom de novo**
8. **ThermoMPNN ΔΔG → Tm 精准预测**
9. **SaProt 结构感知 PLM 评分**

### 必交材料

10. **设计思路 PDF 文档**
11. **GitHub 开源仓库 + README**

---

## 九、关键文件索引

### 提交文件
- `work/round5/submission_round5_hybrid.csv` — 当前推荐提交
- `work/round5/final_6_round5_hybrid.json` — Top-6 详细数据

### 文献
- `work/round5/literature/round5_literature_kb.md` — 完整文献库 (33+ 篇)
- `work/round5/literature/PAPERS_TO_DOWNLOAD.md` — 受订阅限制清单
- `work/round5/文献/` — 4 篇已下载 PDF

### 代码
- `work/round5/01_ligandmpnn_design.py` — LigandMPNN 初版
- `work/round5/03_ligandmpnn_v2.py` — LigandMPNN 修正版
- `work/round5/08_lmpnn_expanded.py` — 多温度扩展
- `work/round5/09b_evolvepro_xgb.py` — EVOLVEpro 训练
- `work/round5/13_evolvepro_score.py` — EVOLVEpro 评分
- `work/round5/17_corrected_selection.py` — 修正版选择

### 数据
- `work/round5/evolvepro_xgb.model` — XGBoost 模型
- `work/round5/evolvepro_scored.json` — 104 条候选含 EVOLVEpro 评分
- `work/round5/esmfold_lmpnn_v2.json` — LigandMPNN v2 ESMFold 结果
- `work/round5/esmfold_lmpnn_expanded.json` — 扩展 LigandMPNN ESMFold 结果
- `work/round5/lmpnn_expanded/` — 190 条 LigandMPNN 原始输出

### 工具
- `work/round5/LigandMPNN/` — LigandMPNN 仓库 + 模型权重
- `work/round5/GeoEvoBuilder/` — GeoEvoBuilder 仓库（缺权重）
- `work/round5/pdbs/` — PDB 结构 (2B3P, 2WUR, 1GFL)

---

## 十、环境信息

```
Python: 3.14.0
PyTorch: 2.11.0+cu128
GPU: NVIDIA GeForce RTX 5080 Laptop GPU (16GB, sm_120)
CUDA: 12.8
OS: Windows 11

关键依赖:
- transformers 4.48.1 (ESMFold)
- xgboost 3.3.0 GPU
- ProDy 2.6.1 (LigandMPNN)
- biopython (LigandMPNN)
- ESM-2 650M: facebook/esm2_t33_650M_UR50D (hf-mirror 下载)

已知问题:
- np.int 已在 LigandMPNN/openfold 中修复
- ESMFold 模型需 local_files_only=True (SSL 问题)
- 中文路径需用 pathlib.Path(r"...") 或英文临时目录
- HuggingFace 需设置 HF_ENDPOINT=https://hf-mirror.com
```

---

## 十一、总结

Round 5 的核心价值不在于生成了更好的候选（没有），而在于：

1. **搭建了 LigandMPNN + EVOLVEpro 工具链**（下一个 Agent 可直接用）
2. **完成了 33+ 篇前沿文献调研**（明确了 GeoEvoBuilder/ESM3/htFuncLib 是下一步关键）
3. **确认了 EVOLVEpro/XGBoost 在 OOD 上仍然失效**（val R² 0.5 但预测偏低）
4. **确认了 mBaoJin 在 ESMFold 中不可用**（pTM 0.36，需要 ColabFold/Boltz-1 验证）

**当前最佳候选仍是 Round 4 的 MPNN_T01_014** (pLDDT 68.3, pTM 0.765)。

**下一个 Agent 应优先做**: GeoEvoBuilder（需下载权重）+ ColabFold 交叉验证 + htFuncLib 16K 数据 fine-tune。

---

*文档作者: Trae AI Agent (Claude)*
*最后更新: 2026-06-22 18:00*
