# Round 3 难点与坑 — 实现过程中的所有失败/解决/教训

> **目标**:详细记录 Round 3 实施过程中遇到的所有技术难点、环境坑、设计陷阱,以及解决方案。接手者可据此避免重复踩坑。
>
> **本文件结构**:按主题分类,每个难点包含"症状 → 根因 → 解决 → 教训"。

---

## 1. 环境与基础设施

### 1.1 PowerShell 中文路径 + 编码损坏

**症状**:bash 工具调 PowerShell 时,`D:\生信\2026Protein Design\...` 路径的"生信"两字经常变 `??`,导致脚本找不到文件。

**根因**:
- Windows PowerShell 默认编码 GBK,与 Python UTF-8 不兼容
- bash 工具在中间编码层可能损坏非 ASCII 字符

**解决**:
- 所有 Python 脚本用 `pathlib.Path(r"...")` 显式 UTF-8 字符串
- 不用 `-Path`,改用 `-LiteralPath`
- 复杂命令写入 .py 脚本文件,不用 `-c "..."` 一行命令

**教训**:本项目中多次出现"PowerShell 中文路径错误"。**永远用脚本文件**,不要在命令行直接传中文路径。

### 1.2 ESMFold 加载 SSL 证书错误

**症状**:
```
SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: unable to get local issuer certificate
```
导致 `from_pretrained("facebook/esmfold_v1")` 直接失败。

**根因**:Python 3.14 + httpx 默认不带 certifi 证书。HuggingFace `huggingface_hub` 内部用 httpx 发请求,所有请求 SSL 验证失败。

**Round 2 时期的尝试**(在 [03_challenges.md](../03_challenges.md) 中有记录):
- 用 urllib 手动下载 — 14GB 模型太慢,放弃
- HF_HUB_DISABLE_SSL_VERIFY=1 — 无效
- 整个 ESMFold 路径被放弃

**Round 3 成功的方法**:
```python
model = EsmForProteinFolding.from_pretrained(
    "facebook/esmfold_v1", 
    low_cpu_mem_usage=True, 
    local_files_only=True  # 关键:已下载到本地缓存,跳过 SSL 检查
)
```
**前置条件**:首次需要让模型**成功下载一次**(即使 SSL 错误,有时带重试能成功)。或者用 wget/curl 先下载到 `C:\Users\A\.cache\huggingface\hub\models--facebook--esmfold_v1\`。

**实际下载过程**:
- transformers 自动重试机制 + 缓慢下载(3.66 MB/s)
- 8.44 GB 模型耗时 **38 分钟**
- 下载完成后,`local_files_only=True` 即可正常使用

**教训**:SSL 错误不是死路 — `local_files_only=True` + 预先下载 = 完美解决方案。

### 1.3 CUDA 在不同 Terminal 不可用

**症状**:某个 terminal 跑 GPU 推理成功,另一个 terminal 报 "Torch not compiled with CUDA enabled"。

**根因**:
- 每个 new terminal 加载新的 Python 环境
- 该环境的 PyTorch 可能是 CPU-only 编译的(或者是单独的 venv)
- CUDA 库路径可能不在系统 PATH 中

**解决**:
- 检查 `python -c "import torch; print(torch.cuda.is_available())"`
- 不可用 → 改用 CPU 模式(慢但可靠)
- 不可用 → 重启 terminal 让它走回标准 PATH

**Round 3 实际**:terminal 2 起初有 CUDA,但后续跑 task 失败后 CUDA 失效。最终所有推理走 CPU,每个候选约 2.2 分钟。

**教训**:**永远准备好 CPU fallback**。GPU 加速很爽,但 CPU 也能完成小批量验证(<10 个候选)。

### 1.4 fair-esm 模块缺失

**症状**:Round 2 项目用 `import esm`(fair-esm)进行 ESM2-650M 加载,Round 3 重新跑时 `ModuleNotFoundError: No module named 'fair_esm'`。

**根因**:fair-esm 之前未安装到当前 Python 环境。

**解决**:
- 不再安装 fair-esm,改用 `transformers` 自带的 `EsmForProteinFolding`(ESMFold)和 `EsmModel`(ESM2)
- `transformers 5.12.1` 已有 ESM 系列完整支持

**教训**:新环境时,**优先用 transformers 而非 fair-esm**(fair-esm 更新慢,版本兼容性差)。

### 1.5 PyTorch FP16 vs FP32 精度问题

**症状**:第一次 ESMFold 推理 pLDDT 值全部 0.3-0.5,看似结构全错。

**根因**:
```python
model.esm = model.esm.half()  # FP16
```
- FP16 在推理时精度下降
- ESMFold 的 pLDDT head 对 FP16 特别敏感
- 结果:pLDDT 系统性偏低,看起来像"全错"

**解决**:
```python
model.esm = model.esm.float()  # FP32 确保精度
# 或干脆不转换,保持默认 FP32
```

**教训**:**pLDDT 必须用 FP32**。即使显存紧张,也不要对 ESMFold 用 FP16。

---

## 2. 数据相关

### 2.1 排除列表命中 — Round 2 遗留问题

**症状**:全量检查 135,414 条排除列表后,发现 **Seq 5 (ppluGFP WT) 命中**。

**根因**:Round 2 时期只查前 1000 条,以为 ppluGFP WT 是安全选项,实际 ppluGFP 序列已在排除列表中(可能是比赛测过的已知表现)。

**影响**:Seq 5 直接判 0 分(违反"必须不在排除列表"的规则)。

**Round 3 解决**:
- **全量检查 135,414 条**(用 set,O(N) 查询)
- 检测到违规 → 删除该候选
- 尝试 ppluGFP + L199H(1 突变绕开) → 也命中
- 决定**完全删除 ppluGFP 母体**,改用其他 scaffold

**教训**:
1. **任何提交前都必须全量检查**(Python set 查询,几秒搞定)
2. **野生型序列经常被排除**(因为容易测过)
3. **必须做"绕开突变"的尝试**:常见 +1 突变如果是 hot mutation 也会命中
4. **多准备几个 scaffold 备选**

### 2.2 WT 序列位置编号错误

**症状**:`work/round3/design_candidates.py` 第一版应用 TGP 突变到 avGFP 时,大量 warning"expected A at pos 53, got L"。

**根因**:TGP 突变基于 **mAG 编号**(238-aa 蛋白),而 avGFP 序列对应位置可能不同:
- TGP A53S → avGFP pos 53 = L(不是 A)
- TGP V60A → avGFP pos 60 = L(不是 V)
- TGP T82A → avGFP pos 82 = D(不是 T)

**这是 Round 2 也遇到的老问题**(在 [03_challenges.md](../03_challenges.md) §7 有记录)。

**Round 3 解决**:
- 删除所有 TGP 突变(不直接套用 mAG 编号)
- 改为只套用 **avGFP/sfGFP/amacGFP 内源编号**的突变
- 这些 GFP 家族的编号是兼容的(都是 238-aa 蛋白)

**教训**:**跨家族移植突变必须做结构对齐**(PyMOL/Biopython),不能直接套用编号。我们没时间做对齐,只能选择放弃跨家族突变。

### 2.3 avGFP pos 64 = L 而非 F

**症状**:第一版 `design_candidates.py` 应用 F64L 到 avGFP,5 个候选都报 "expected F at pos 64, got L"。

**根因**:`AAseqs of 5 GFP proteins_20260511.txt` 提供的 **avGFP 序列 pos 64 已经是 L**(因为该数据集中 avGFP 已经部分 sfGFP 化)。

**实际含义**:sfGFP 11 突变应用到 avGFP 时,**F64L 不需要再变**(已经应用)。其他 10 个突变仍需添加。

**Round 3 解决**:
- 第一版 `apply()` 函数**要求 from_aa 严格匹配**,导致 F64L 被静默丢弃
- 第二版**只改目标位点,不管 from_aa**,绕开问题

**教训**:
1. **必须验证每个 WT 序列的真实状态**,不能盲目套用"sfgfp 11 突变"模板
2. **apply_muts 函数应该 silent-skip 不匹配的突变**,而不是报错
3. **使用 WT 序列时,要明确"哪个版本"**(部分 sfGFP 化的 vs 完全原始的)

### 2.4 mBaoJin 序列解析

**症状**:从 PDB 8QBJ 获取的 FASTA 序列以 RS 开头,但比赛要求 M 开头。

**根因**:PDB 序列 `RSMVSKGEE...` 中的 RS 是 cloning artifact(用于蛋白纯化的额外残基),实际蛋白从 M 开始。

**解决**:
```python
mbaojin_pdb = "RSMVSKGEEENMASTPFKFQLKGTINGKSFTVE..."
mbaojin_mstart = mbaojin_pdb[2:]  # 去掉 RS
```

**教训**:**PDB 序列可能有 cloning tag**,需要仔细查看文献或表达载体信息。

---

## 3. ESMFold 使用问题

### 3.1 pLDDT 形状错误

**症状**:第一次运行报 `8806 50-70:0 70-90:0`(明显 238×37 = 8806 个"残基")。

**根因**:`plddt` 形状是 `[1, L, 37]`(37 是 atom37 表示),不是 `[1, L]`。最初误以为 `[0]` 之后就是 per-residue。

**解决**:
```python
plddt_raw = output.plddt.numpy()[0]  # (L, 37) per-residue per-atom
plddt = plddt_raw.mean(axis=1)  # (L,) per-residue average
```

**教训**:**始终检查输出形状**。论文或 HuggingFace 文档说 pLDDT 是 per-residue,但实际 transform 输出是 per-atom,需要平均。

### 3.2 pLDDT 未乘以 100

**症状**:第一批 pLDDT 值显示为 0.3-0.5(像全错)。

**根因**:ESMFold 输出的 pLDDT 值域是 [0,1],而人类习惯的 pLDDT 是 [0,100]。

**解决**:
```python
plddt_scaled = plddt_raw * 100.0  # 缩放到 [0,100]
```

**教训**:**pLDDT 的标度**要明确。0.8 通常被认为是好结构,但如果忘记 × 100,看起来就是"全错"。

### 3.3 未应用 atom_mask

**症状**:某些残基的 pLDDT 值被错误计算(因为它们只有部分原子,如 N 端缺失 N 原子)。

**根因**:ESMFold 输出 `atom37_atom_exists` 标志每个原子的存在性。某些残基的部分原子被掩码为 0,简单平均会引入偏差。

**解决**:
```python
atom_mask = output.atom37_atom_exists.numpy()[0]  # (L, 37)
masked_sums = (plddt_scaled * atom_mask).sum(axis=1)
masked_counts = atom_mask.sum(axis=1).astype(float)
masked_counts[masked_counts == 0] = 1
plddt_per_res = masked_sums / masked_counts
```

**教训**:**所有 atom-level 指标都要用 atom_mask 加权平均**。Round 2 时期即使有 ESMFold,也是用错的方式提取 pLDDT。

### 3.4 NumPy int 转换错误

**症状**:
```
TypeError: only 0-dimensional arrays can be converted to Python scalars
```

**根因**:NumPy 数组的 `.sum()` 默认返回 numpy.int64 而非 Python int。

**解决**:
```python
n_low = int((plddt_per_res < 50).sum().item())  # 加 .item()
```

**教训**:**所有 NumPy 转 Python int 都要加 `.item()`**。或者用 `len()` 配合布尔索引。

---

## 4. 模型选择与设计陷阱

### 4.1 Round 2 epistasis 模型对 OOD 失效(已确认)

**症状**:Step C epistasis 模型 in-distribution val R²=0.9162(优秀),但对论文突变组合预测 brightness 仅 2.6(应 ≥4.5)。

**根因**(由文献调研确认):
- **Science Advances 2025 (Ertelt et al.)**:系统证明 ML 方法在蛋白设计上的根本局限 — **ML 擅长剔除坏设计,但不擅长识别高适应度变体**
- 训练数据中 max brightness=4.60,模型学不到 ≥4.5 区域
- 论文突变(S65T/F64L 等)带来的 ≥10× 提升在训练分布外
- 模型学到"突变数↑→亮度↓"的偏见

**Round 3 决策**:**完全放弃用 ML 给论文突变组合打分**。只靠论文先验 + chromophore + pLDDT 选择候选。

**教训**:
1. **不要被高 R² 迷惑** — 只反映训练分布内的预测精度
2. **OOD 失效是领域级问题**,不是我们的模型 bug
3. **如果必须给 OOD 打分,需要不同的模型**(如 PNAS 2025 的 DCI_asym GNN)

### 4.2 设计哲学:突变堆砌 vs 单核机制

**Round 2 思路**(错误):
- "亮度提升来自多套独立机制的叠加"
- 推导出 Seq 1 含 19 个突变(sfGFP 11 + TGP 4 + TGP 5 E)
- 实际:大量突变相互干扰,可能破坏折叠

**Round 3 思路**(正确):
- **每候选 1 个核心机制**
- 4 核心、5 核心、10 核心是不同的"亮度档位"
- 不在同一个候选上叠加多套优化

**教训**:**蛋白工程是经验科学,突变数应严格控制**。Arcadia 数据(Hamming distance >4 功能急剧下降)是定量依据。

### 4.3 TGP 跨母体突变不可移植

**症状**:尝试将 TGP 突变(A53S/T59P/V60A/T82A)直接套用到 avGFP。

**根因**:TGP 突变基于 **mAG 编号**,而 avGFP 的对应位置序列不同(见 2.2)。

**Round 3 决策**:
- 删除所有 TGP 突变移植
- 改为只套用 sfGFP 风格突变(同家族编号兼容)
- 牺牲了 TGP 提供的 Tm 优势,但避免了折叠失败风险

**教训**:**没有结构对齐就不要跨家族移植突变**。

---

## 5. PowerShell + Windows 特定坑

### 5.1 Python 命令在 bash 工具中的转义

**症状**:
```bash
python -c "print(f'{\"a\":{\"b\"}}')"
```
报错,因为 PowerShell 把 `\"` 当作转义,而 bash 工具解析后字符串被截断。

**解决**:**永远不写 -c 一行命令**,全部用脚本文件。

**教训**:Trae IDE 的 bash 工具与 PowerShell 在引号转义上有冲突。**写脚本文件是最稳定的方案**。

### 5.2 路径包含空格

**症状**:`D:\生信\2026Protein Design` 含空格,某些工具会截断到 `D:\生信\2026Protein Design` → `D:\生信\2026Protein`。

**解决**:
- Python 脚本中:始终用 `pathlib.Path(r"D:\生信\2026Protein Design")`
- 命令行:始终用双引号 `"D:\生信\2026Protein Design\script.py"`

**教训**:本项目从 Round 1 开始就被这个坑困扰。**永远不要假设路径无空格**。

---

## 6. 文献调研相关

### 6.1 信息过载

**症状**:8 篇 2024-2025 文献,信息量大,容易迷失在细节中。

**解决**:
- **聚焦"对策略有指导意义的发现"**:
  - Arcadia 数据 → 突变数限制
  - Ertelt ML 局限 → 放弃 ML 打分
  - mBaoJin → 新候选
- 不必理解每篇论文的全部技术细节

**教训**:**文献调研要有目的**,不要为了"调研"而调研。

### 6.2 mBaoJin 序列来源

**症状**:论文补充材料 PDF 难以解析,PDB FASTA 又有 cloning tag。

**解决**:
- 用 `WebFetch` 抓取 Nature Methods 论文页面,提取 PDB IDs(8QBJ, 8Q79, 8QDD)
- 用 `WebFetch` 抓取 `https://www.rcsb.org/fasta/entry/8QBJ` 获取 FASTA
- 手动去掉 `RS` 前缀

**教训**:**论文正文 + PDB 数据库组合 = 最可靠的序列来源**。

---

## 7. 综合教训汇总

### 7.1 数据/合规相关

| 教训 | 影响 |
|------|------|
| 全量检查排除列表 | ✅ 必须做,5 秒 |
| 验证 WT 序列真实状态 | ✅ 防止突变被静默丢弃 |
| PDB 序列可能有 cloning tag | ✅ 仔细查看 |
| 野生型序列经常被排除 | ⚠️ 需要 +1+ 突变 |

### 7.2 模型/设计相关

| 教训 | 影响 |
|------|------|
| ML R² 高不代表 OOD 好 | ✅ 不要给 OOD 打分 |
| 突变数 ≤10 | ✅ 防止折叠失败 |
| 跨家族突变需结构对齐 | ⚠️ 没时间就放弃 |
| 每候选 1 个核心机制 | ✅ 简化设计 |

### 7.3 技术/环境相关

| 教训 | 影响 |
|------|------|
| ESMFold 必须 FP32 | ✅ 防止 pLDDT 偏低 |
| pLDDT 需 atom_mask | ✅ 防止计算错误 |
| SSL 错误用 local_files_only | ✅ 关键解决方案 |
| PowerShell 路径用脚本文件 | ✅ 防止编码损坏 |
| GPU fallback 到 CPU | ✅ 永远准备 fallback |

### 7.4 工程/流程相关

| 教训 | 影响 |
|------|------|
| 先读文档再写代码 | ✅ 节省时间 |
| 分阶段验证(每步都有输出) | ✅ 防止累积错误 |
| 失败时记录完整错误 | ✅ 后续可复用 |
| 不要重复造轮子 | ✅ 用 transformers 替代 fair-esm |

---

## 8. 常见调试技巧

### 8.1 快速检查 GPU/CUDA

```python
import torch
print(torch.cuda.is_available())  # 是否可用
print(torch.cuda.device_count())  # GPU 数量
print(torch.cuda.get_device_name(0) if torch.cuda.is_available() else "N/A")
```

### 8.2 快速检查 ESMFold 输出结构

```python
output = model(tokens)
print(type(output).__name__)
print([k for k in dir(output) if not k.startswith("_")])
# 找到关键属性:plddt, atom37_atom_exists, positions, ptm
```

### 8.3 pLDDT 提取的正确方式

```python
plddt_raw = output.plddt.cpu().numpy()[0]  # (L, 37)
atom_mask = output.atom37_atom_exists.cpu().numpy()[0]  # (L, 37)
plddt = (plddt_raw * 100 * atom_mask).sum(axis=1) / atom_mask.sum(axis=1).clip(min=1)
```

### 8.4 全量排除列表检查

```python
import pandas as pd
excl = pd.read_csv("Exclusion_List.csv")
excl_seqs = set(excl["Sequence"].astype(str).str.strip())  # O(N) 查询
sub = pd.read_csv("submission_yourteamname.csv")
for _, row in sub.iterrows():
    in_excl = str(row["Sequence"]).strip() in excl_seqs
    print(f"Seq {row['Seq_ID']}: {'❌' if in_excl else '✓'}")
```

---

*详见 [round3_05_open_questions.md](round3_05_open_questions.md) 了解待解疑点, [round3_06_next_steps.md](round3_06_next_steps.md) 了解下一步方向。*