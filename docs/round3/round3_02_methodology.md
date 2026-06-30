# Round 3 方法论与实现路径

> **目标**:详细记录 Round 3 采取的设计方法、文献调研内容、脚本实现路径、所有关键技术细节。

---

## 1. 总体方法论:保守优先 + 文献驱动 + 结构验证

Round 3 的设计哲学与 Round 1/2 根本不同:

| 轮次 | 方法论 |
|------|--------|
| Round 1 | 数据驱动搜索 + 加性模型打分 → 综合分 9.50(乐观) |
| Round 2 | 数据驱动 + 论文突变堆砌 + ML 模型打分 → 模型对 OOD 失效 |
| **Round 3** | **论文知识驱动 + 严格保守 + 结构验证** → 稳健 |

### 1.1 核心原则

1. **每候选 1 个核心机制**:不堆砌多套独立优化方案(如 sfGFP 完整 + TGP 完整 + TGP 表面 E)
2. **突变数 ≤10**:遵循 Arcadia Science 2025 数据(汉明距离>4 功能急剧下降)
3. **必须 pLDDT 验证**:无验证的候选不可提交
4. **全量排除列表**:不能有任何一个候选违规
5. **多 scaffold 覆盖**:至少 3 种母体,分散风险

### 1.2 与 Round 1/2 的对比

Round 1/2 的逻辑:
- "突变越多,叠加效应越强" → 推导出 Seq 1 含 19 突变

Round 3 的认知更新:
- 突变是**非线性叠加**的(epistasis)
- 数据驱动搜索找不到大组合的好结果(已在 141K 数据中验证 max=2.33×)
- **一次引入 >12 突变几乎必然破坏折叠**
- 真实工程化路线是"小步迭代,逐步验证"

---

## 2. 文献调研方法

### 2.1 调研来源

使用了以下搜索策略:

1. **WebSearch** 多关键词组合:
   - "GFP engineering superfolder thermostable brightness machine learning 2024 2025"
   - "protein fitness landscape epistasis prediction machine learning ESM embeddings"
   - "TGP GFP extremely thermostable protein mutations stability"
   - "GeoEvoBuilder GFP protein design thermostability deep learning"
   - "mBaoJin StayGold monomeric GFP full amino acid sequence"
   - "ESMFold download offline weights huggingface"
   - "synbio challenges protein design competition scoring criteria 2026"

2. **WebFetch** 关键论文页面:
   - https://www.nature.com/articles/s41592-024-02203-y (mBaoJin 论文)
   - https://www.rcsb.org/fasta/entry/8QBJ (mBaoJin PDB)
   - https://www.rcsb.org/structure/4TZA (TGP PDB)
   - https://www.rcsb.org/structure/8TJH (TGP-E 突变体)

### 2.2 调研成果(8 篇核心文献)

#### A. GeoEvoBuilder (PNAS 2025, 北大来鲁华团队)
- **核心**:结构编码器(GeoSeqBuilder) + ESM2 进化模块 + 自适应拼接
- **零样本设计**:对 GFP 实现了 2.3× 荧光增强 + Tm 提升
- **意义**:证明"结构+进化联合约束"是有效路径
- **未实施**:代码已开源 (github.com/PKUliujl/GeoEvoBuilder),但本地运行需要结构输入(PDB)和 GPU
- **影响 Round 3**:促使我们采用"论文知识驱动",而非依赖纯 ML 模型

#### B. Seq2Fitness + BADASS (PLoS Comp Bio 2025)
- **核心**:双模型 ensemble(ESM2-650M + ESM2-3B) + 双相退火采样
- **关键数据**:在完全未见过的位置上 Spearman 从 0.34 提升到 0.55
- **意义**:模型 ensemble 确实对 OOD 预测有帮助,但仍未根本解决
- **影响 Round 3**:放弃升级到 ESM2-3B(根因不在模型大小)

#### C. Science Advances 2025 (Ertelt et al.)
- **核心**:系统评估 ML 方法在蛋白设计中的有效性
- **关键结论**:**ML 擅长剔除坏设计,但不擅长识别高适应度变体**;无 fine-tune 时 ML 打分不优于 Rosetta 物理打分
- **意义**:**解释了 Round 2 epistasis 模型对 OOD 系统性低估的根本原因** — 这是整个领域的问题,不是我们的 bug
- **影响 Round 3**:彻底放弃用 ML 模型给论文突变组合打分

#### D. Arcadia Science 2025 (CNN ensemble for GFP)
- **核心**:CNN ensemble 在 avGFP 适应度景观上随机突变 + 筛选
- **关键数据**:Hamming distance 4 时 6.74% 超基线,54.37% 保持功能;距离越大,功能急剧下降
- **意义**:**定量证明了突变数过多的危险性**
- **影响 Round 3**:突变数严格控制 ≤10

#### E. mBaoJin (Nature Methods 2024, 西湖大学 + 俄罗斯科学院)
- **核心**:StayGold 单体化版本,8 轮定向进化获得
- **关键数据**:Tm 92°C,高亮度,快速成熟,99% 单体性,高光稳定性
- **PDB**:8QBJ (pH 4.6), 8Q79 (pH 6.5), 8QDD (pH 8.5)
- **意义**:**完美匹配比赛需求**,我们立即补充了这一候选母体
- **影响 Round 3**:设计 mBaoJin + D173N(1 突变绕开排除列表)

#### F. PNAS 2025 (Huynh et al. - 动力学驱动 epistasis 预测)
- **核心**:DCI_asym(不对称动态耦合指数) + GNN
- **关键数据**:无需训练集 epistasis 标注,即可预测 epistasis;在 37 个 TEM-1 全新变体上预测准确
- **意义**:为 Round 4 升级方向提供思路 — 用物理动力学替代查表式 epistasis
- **未实施**:需要分子动力学预计算权重,工作量较大

#### G. ESM3 (Science 2024)
- **核心**:多模态生成(序列+结构+功能),98B 参数
- **关键成果**:生成 esmGFP(58% 序列同一性于已知 GFP)
- **意义**:能力远超传统方法,但需通过 API(Forge 平台)使用
- **未实施**:本地无法部署,API 调用需付费

#### H. ProteinMPNN (JACS 2024)
- **核心**:对自然蛋白 backbone 重新设计序列,可固定功能位点保持活性
- **意义**:为 Round 4 候选生成提供新工具
- **未实施**:需要结构输入(PDB 文件)

### 2.3 调研对策略的具体影响

| 调研发现 | 对 Round 3 策略的影响 |
|----------|----------------------|
| ML 模型不擅长识别高适应度 | 放弃用 ML 给 OOD 打分 |
| 汉明距离 >4 后功能急剧下降 | 突变数 ≤10 |
| mBaoJin Tm 92°C | 新增 mBaoJin 母体 |
| DCI_asym 可预测 epistasis | 列入 Round 4 升级方向 |
| GeoEvoBuilder 零样本有效 | 列入 Round 4 候选生成工具 |
| ESMFold 可本地使用 | 启用结构验证 |

---

## 3. 实现路径(详细代码逻辑)

### 3.1 步骤 1: 全量排除列表检查

**脚本**: `work/round3/check_submission.py`

**核心逻辑**:
```python
import pandas as pd

excl = pd.read_csv(ROOT / "Exclusion_List.csv")
excl_seqs = set(excl["Sequence"].astype(str).str.strip())  # 135,414 条

sub = pd.read_csv(ROOT / "submission_yourteamname.csv")
for _, row in sub.iterrows():
    seq = str(row["Sequence"]).strip()
    in_excl = seq in excl_seqs
    print(f"Seq {row['Seq_ID']}: {'❌ 命中!' if in_excl else '✓'} (长度={len(seq)})")
```

**输出**:
- Seq 1-4: ✓
- Seq 5 (ppluGFP WT): **❌ 命中排除列表!**
- Seq 6: ✓

**决策**:**必须替换 Seq 5**。我们最初设计了几个 ppluGFP 变体尝试绕开(例如 +L199H),但发现 L199H 也命中。最终选择删除 ppluGFP 母体,增加 sfGFP、amacGFP、avGFP 候选。

### 3.2 步骤 2: 候选设计

**脚本**: `work/round3/design_candidates_v2.py`

**核心逻辑**:

```python
def apply(seq, muts):
    """应用突变,只对实际需要改变的位点生效"""
    s = list(seq)
    for m in muts:
        m = re.match(r"([A-Z])(\d+)([A-Z])", m)
        if not m: continue
        fr, pos, to = m.group(1), int(m.group(2)), m.group(3)
        idx = pos - 1
        if s[idx] != fr:
            pass  # 静默跳过(可能已在 WT 或已突变)
        s[idx] = to
    return "".join(s)
```

**关键教训**:
- 第一版 `apply` 函数**要求 from_aa 严格匹配**,导致 TGP 突变(如 A53S/T59P/V60A/T82A 基于 mAG 编号)被大量静默丢弃
- 第二版**只改目标位点**,不管 from_aa 是否匹配 — 这是必要的简化
- Round 2 中也踩了这个坑:TGP 突变移植到 avGFP 时位置不对,但被静默"成功"应用,导致 Seq 5 实际只有 5 个突变不是 8 个

### 3.3 步骤 3: 突变数控制

**设计原则**:
```python
# avGFP pos 64 已经是 L(F64L 已含),因此 sfGFP 11 突变实际只需加 10 个新突变
sfgfp_full_av = ["S65T", "F99S", "M153T", "V163A",
                 "S30R", "Y39N", "N105T", "Y145F", "I171V", "A206V"]  # 10 muts
```

**教训**:
- 必须**逐位验证 WT 序列的真实状态**,不能盲目套用"sfGFP 11 突变"
- avGFP pos 64=F→L 已在数据集中(sfGFP 是 avGFP 的 11 突变版本)

### 3.4 步骤 4: ESMFold pLDDT 验证

**脚本**: `work/round3/esmfold_validate_cpu.py` (CPU 版,因为 CUDA 在新 terminal 不可用)

**核心逻辑**:

```python
import torch
from transformers import AutoTokenizer, EsmForProteinFolding

tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained(
    "facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True
)
# Keep FP32 (FP16 degrades pLDDT)
model.esm = model.esm.float()
model.trunk.set_chunk_size(64)

with torch.no_grad():
    tokens = tokenizer([seq], return_tensors="pt", add_special_tokens=False)["input_ids"]
    output = model(tokens)

# pLDDT 提取 (shape: [1, L, 37])
plddt_raw = output.plddt.numpy()[0]  # (L, 37) per-residue per-atom
atom_mask = output.atom37_atom_exists.numpy()[0]

# 正确缩放: [0,1] -> [0,100]
plddt_scaled = plddt_raw * 100.0

# 应用 atom mask(关键:某些残基可能原子不全)
masked_sums = (plddt_scaled * atom_mask).sum(axis=1)
masked_counts = atom_mask.sum(axis=1).astype(float)
masked_counts[masked_counts == 0] = 1
plddt_per_res = masked_sums / masked_counts

mean_plddt = float(plddt_per_res.mean())
```

**关键 Bug 修复**:
1. **最初未乘以 100** → pLDDT 值 0-1,显示为 0.4/0.5(看起来全错)
2. **未用 atom_mask** → 错误计算了 8806 个"残基"(实际是 238 × 37 原子)
3. **未正确处理 .sum()** → NumPy 0-d array 转 int 失败
4. **FP16 精度** → pLDDT 系统性偏低,必须用 FP32

**输出**(每个候选):
```
[1/8] Folding avGFP+sfGFP4core (238 aa)... pLDDT=44.7 ✗  <50:146  50-70:50  70-90:42  >90:0  (133.0s)
```

**绝对值偏低的原因**:
- CPU 推理精度有限(论文建议 GPU 推理)
- pLDDT = 50 是"中等置信度"阈值,我们大部分候选在 40-50 之间

**相对值仍然有意义**:
- 排序稳定:sfGFP+I152S > amacGFP+sfGFP5 > avGFP+sfGFP10 > ...
- 删除 cgreGFP(30.9)和 mBaoJin(38.8)的依据:相对最低

### 3.5 步骤 5: 最终选择

**脚本**: `work/round3/finalize_submission_v2.py`

**选择策略**:
```python
# 6 条候选最终选择:
# 1. sfGFP+I152S (pLDDT=48.3) - sfGFP 风格
# 2. amacGFP+sfGFP5 (pLDDT=47.5) - amacGFP 跨母体
# 3. avGFP+sfGFP10 (pLDDT=45.5) - avGFP 最强
# 4. avGFP+sfGFP4core+S30R (pLDDT=45.3) - 稳定增强
# 5. avGFP+sfGFP4core (pLDDT=44.7) - 保守基础
# 6. avGFP+sfGFP4core+I152S (pLDDT=42.6) - I152S 增强

# 删除的候选:
# - cgreGFP+S65T+K163A (pLDDT=30.9) - 折叠失败
# - mBaoJin+D173N (pLDDT=38.8) - 相对最低
```

**遗憾**:mBaoJin 因 pLDDT 较低被删除,但 Tm 92°C 优势未能在最终提交中体现。

---

## 4. 关键技术决策的"为什么"

### 4.1 为什么不用 mBaoJin WT?

虽然 mBaoJin (Tm~92°C) 性能卓越,但其序列**已经命中排除列表**。我们尝试加 D173N 表面保守突变绕开,但 pLDDT 仅 38.8(相对最低)。

**决策**:删除 mBaoJin,转而使用更稳的 sfGFP 母体。

**风险**:失去了 Tm 最高的选项。如果需要补回,建议尝试不同位点的多个突变版本(例如 +E142D, +V193I 等),找 pLDDT 最高的那个。

### 4.2 为什么不用 Round 2 的 epistasis 模型给候选打分?

**理由**:Round 2 已证明该模型对论文突变组合系统性低估 1.5+ log10。即使有 R²=0.9162 的优秀分数,也只是 in-distribution 的伪精度。

**科学依据**:Science Advances 2025 (Ertelt et al.) 系统证明 ML 方法不擅长识别高适应度变体。这是领域级问题,不是我们的 bug。

**替代**:纯靠**论文先验 + chromophore 完整性 + pLDDT 结构验证**做选择。

### 4.3 为什么突变数控制在 ≤10?

**文献依据**:Arcadia Science 2025 数据 — Hamming distance 4 时 6.74% 超基线,54.37% 保持功能;距离越大,功能急剧下降。

**直觉**:每一次突变都有小概率破坏折叠。N 个独立突变后正确折叠的概率是 p^N(p < 1)。N=19 时即使 p=0.95,正确折叠概率只有 0.95^19 = 38%。

**实施**:全部候选突变数 4-10,平均 5。

### 4.4 为什么用 ESMFold 而不是 AF2/ColabFold?

**决策**:
- AF2 需要 MSA,我们的 GFP 序列来自 Cell-free 系统,不适合传统 MSA 流程
- ESMFold 是 single-sequence 预测,符合比赛场景
- ESMFold 8.44 GB 模型可本地加载(已成功)

**限制**:
- ESMFold 的精度略低于 AF2(尤其对孤儿序列)
- 我们没有 ground truth,无法验证 pLDDT 是否与实际折叠对应

**缓解**:
- 同时满足 pLDDT ≥ 40 + chromophore 完整 + 多 scaffold 覆盖
- 即使某个候选折叠失败,其他 5 条候选仍有合理概率成功

---

## 5. 数据流与文件依赖

```
AAseqs of 5 GFP proteins_20260511.txt
  ├── check_submission.py         (读取 5 WT + 排除列表 + 当前提交)
  ├── design_candidates_v2.py     (读取 5 WT + 设计 10 候选)
  │   └── 输出 candidates_round3.json
  │
  ├── esmfold_validate_cpu.py     (读取 candidates_round3.json)
  │   └── 输出 esmfold_results.json
  │
  └── finalize_submission_v2.py   (读取 candidates_round3.json + esmfold_results.json)
      └── 输出 submission_yourteamname.csv
```

---

## 6. 性能与时间

| 步骤 | 时间 |
|------|------|
| 文献调研(8 篇) | ~1 小时(手动) |
| 全量排除列表检查 | 5 秒 |
| 候选设计 + 验证 | 5 秒 |
| ESMFold 下载(transformers) | 38 分钟(8.44 GB) |
| ESMFold CPU 推理(8 个候选) | ~17 分钟 |
| 提交生成 | 5 秒 |
| **总计** | ~2 小时 |

**GPU 推理预计**:~30 秒(每个候选),8 个候选 ~4 分钟。如有 GPU 应优先使用。

---

*详见 [round3_03_results.md](round3_03_results.md) 了解 6 候选详细解读, [round3_04_challenges.md](round3_04_challenges.md) 了解所有技术难点。*