# Round 4 技术参考手册

> **目的**: 提供 Round 4 用到的所有工具、模型、数据集的详细技术规格
> **适用**: 接手 Round 5 复现实验、改进模型、调试问题

---

## 一、硬件环境

| 项目 | 规格 |
|------|------|
| **GPU** | NVIDIA RTX 5080 Laptop |
| **VRAM** | 16 GB |
| **CUDA** | 13.3 (驱动), 12.8 (PyTorch 库) |
| **sm** | 12.0 (Blackwell 架构) |
| **CPU** | (未记录, 但 GPU 推理不依赖) |
| **OS** | Windows |
| **Python** | 3.14.0 |
| **CUDA Toolkit** | 不需要单独装, PyTorch 自带 runtime |
| **NVIDIA 驱动** | 610.47 |

### GPU 推理速度基准

| 任务 | 速度 | 备注 |
|------|------|------|
| ESMFold 1 序列 (238 aa) | 5-6s | GPU FP32, chunk=128 |
| ESMFold 1 序列 (CPU 130s) | 130s | 不推荐 |
| ProteinMPNN 1 序列 (230 aa) | ~5-10s | GPU 推理 |
| ThermoMPNN (未装) | - | 应 ~5s/单点突变 |
| ESM-IF1 (未装) | - | 应 ~10s/序列 |

---

## 二、软件环境

### 2.1 Python 包版本

```
torch                          2.11.0+cu128
torchvision                    0.26.0+cu128
transformers                   4.48.1
esm                           3.2.1.post1 (没用到)
fair-esm                       2.0.0 (没用到)
numpy                          2.4.3
pandas                         (未记录, 应 >= 2.0)
scikit-learn                   (R2 时期, 本轮未直接用)
Pillow                         12.2.0
requests                       2.32+
openpyxl                       (本轮没用到)
matplotlib                     (绘图, 间接)
reportlab                      (PDF 生成, 间接)
```

### 2.2 pip 安装命令

```powershell
# PyTorch GPU (国内推荐用上海交大镜像)
pip install torch torchvision --index-url https://mirror.sjtu.edu.cn/pytorch-wheels/cu128 --trusted-host mirror.sjtu.edu.cn

# Transformers + ESMFold
pip install transformers

# 其他依赖
pip install numpy pandas pillow openpyxl
```

### 2.3 关键库的来源

| 库 | 来源 | 用途 |
|----|------|------|
| ESMFold | HuggingFace `facebook/esmfold_v1` | 3D 结构预测 + pLDDT |
| ProteinMPNN | https://gitcode.com/gh_mirrors/pr/ProteinMPNN | de novo 序列设计 |
| PyTorch | 国内镜像 (上交大/阿里云) | GPU 计算 |
| 5 个 GFP PDB | RCSB PDB | 晶体结构输入 |

---

## 三、ESMFold 详细使用

### 3.1 模型加载

```python
from transformers import AutoTokenizer, EsmForProteinFolding
import torch

tokenizer = AutoTokenizer.from_pretrained(
    "facebook/esmfold_v1", 
    local_files_only=True  # 关键: 已下载, 不重新拉
)
model = EsmForProteinFolding.from_pretrained(
    "facebook/esmfold_v1",
    low_cpu_mem_usage=True,
    local_files_only=True
).cuda()
model.trunk.set_chunk_size(128)  # 长序列必须
model.eval()
```

### 3.2 推理代码

```python
with torch.no_grad():
    tokens = tokenizer(
        [seq], 
        return_tensors="pt", 
        add_special_tokens=False
    )["input_ids"].cuda()
    output = model(tokens)
```

### 3.3 关键输出提取

```python
# pLDDT: 形状 [1, L, 37], 值域 [0, 1]
plddt_raw = output.plddt.cpu().numpy()[0]  # (L, 37)
atom_mask = output.atom37_atom_exists.cpu().numpy()[0]  # (L, 37)

# 缩放到 [0, 100]
plddt_scaled = plddt_raw * 100.0

# 用 atom_mask 加权平均 (避免 N 端缺原子影响)
masked_sums = (plddt_scaled * atom_mask).sum(axis=1)
masked_counts = atom_mask.sum(axis=1).astype(float)
masked_counts[masked_counts == 0] = 1
plddt_per_res = masked_sums / masked_counts  # (L,)

# chromophore 区域 (T65-Y66-G67)
cb_start, cb_end = 58, 73  # 1-based
cb_mean = float(plddt_per_res[cb_start:cb_end].mean())

# pTM (标量, 不是 tensor)
ptm = float(output.ptm.cpu().item())
```

### 3.4 重要参数

| 参数 | 推荐值 | 备注 |
|------|--------|------|
| `chunk_size` | 128 | 减少长序列显存 |
| `FP16 vs FP32` | **FP32** | FP16 严重降低 pLDDT |
| `local_files_only` | True | 已下载避免重连 |
| `add_special_tokens` | False | 不要 BOS/EOS |

### 3.5 常见错误

```python
# ❌ 错误
plddt_raw = output.plddt.numpy()[0]  # shape 错! 应该是 (L, 37) 不是 (L,)
plddt_mean = plddt_raw.mean()  # 这是 atom-level 平均, 不是 residue-level

# ✅ 正确
plddt_per_res = ((plddt_raw * 100) * atom_mask).sum(axis=1) / atom_mask.sum(axis=1)
plddt_mean = plddt_per_res.mean()

# ❌ 错误
ptm = output.ptm[0].cpu().numpy()  # IndexError: ptm 是 0-dim tensor

# ✅ 正确
ptm = output.ptm.cpu().item()
```

---

## 四、ProteinMPNN 详细使用

### 4.1 源码获取

```bash
# 国内推荐
git clone https://gitcode.com/gh_mirrors/pr/ProteinMPNN.git
cd ProteinMPNN

# 切到纯英文路径 (否则会报 FileNotFoundError)
# 将整个 ProteinMPNN 复制到 C:\Temp\mpnn_work2\ProteinMPNN
```

### 4.2 工作流程

```bash
# 1. 解析 PDB
python helper_scripts/parse_multiple_chains.py \
    --input_path ./pdbs_dir \
    --output_path parsed.jsonl

# 2. 指定要设计的链
python helper_scripts/assign_fixed_chains.py \
    --input_path parsed.jsonl \
    --output_path chains.jsonl \
    --chain_list "A"

# 3. 固定关键位
python helper_scripts/make_fixed_positions_dict.py \
    --input_path parsed.jsonl \
    --output_path fixed.jsonl \
    --chain_list "A" \
    --position_list "61 62 63 64 65 66 67 68 69 70 96 145 148 167 203 205 222"

# 4. 跑设计
python protein_mpnn_run.py \
    --jsonl_path parsed.jsonl \
    --chain_id_jsonl chains.jsonl \
    --fixed_positions_jsonl fixed.jsonl \
    --out_folder out_T01 \
    --num_seq_per_target 30 \
    --sampling_temp 0.1 \
    --seed 37 \
    --batch_size 8 \
    --save_score 1 \
    --model_name v_48_020 \
    --path_to_model_weights /path/to/vanilla_model_weights
```

### 4.3 重要参数

| 参数 | 推荐值 | 备注 |
|------|--------|------|
| `sampling_temp` | 0.1 (保守) / 0.3 (多样) | 低温保守, 高温多样 |
| `model_name` | v_48_020 | 平衡点 (0.2Å 噪声) |
| `num_seq_per_target` | 20-30 | 太多会重复 |
| `batch_size` | 8 | 视 GPU 而定 |
| `seed` | 37 | 可重复 |

### 4.4 输出格式

```
>.\2B3P, score=0.8987, global_score=0.9393, fixed_chains=[], designed_chains=['A'], ...
SKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLXXXVQCFSRYPD...
> T=0.1, sample=1, score=0.5742, global_score=0.6347, seq_recovery=0.7477
GKGDELFKGVVPVLVELDGDVNGHKFSVKGEGEGDASKGKLTLKFVCTTGKLPVPWPTLVTTLXXXVQCFAKYPEHMK...
```

**关键**:
- `score` = -log(p), 越**低**越好
- `seq_recovery` = 0-1, 越高越接近原始
- `XXX` = 固定位占位符, 需要用 PDB WT 填回
- 第一条是 PDB 原始序列 (含 X)

### 4.5 X 占位符处理

```python
# 1. 读取 PDB WT
with open(fa_path) as f:
    content = f.read()
entries = content.strip().split(">")
wt_with_x = entries[1].split("\n", 1)[1].strip()  # 第一条是 WT

# 2. 找对齐位置
sfgfp_wt = "MSKGEELFTG..."  # 实际 WT 序列
for offset in range(0, 15):
    # 比对前 30 个非 X 字符
    matches = sum(1 for a, b in zip(wt_with_x[:30], sfgfp_wt[offset:offset+30]) if a == b or a == 'X')
    if matches >= 28:  # 找到对齐
        break

# 3. 填回 X
wt_aligned = sfgfp_wt[offset+1:offset+1+len(wt_with_x)]
filled = "".join(wt_aligned[i] if c == 'X' else c for i, c in enumerate(mpnn_seq))

# 4. 补 N 端
full_seq = sfgfp_wt[:offset+1] + filled + sfgfp_wt[offset+1+len(filled):]
```

### 4.6 常见错误

```
FileNotFoundError: 'C:\\Temp\\...out\\T01//seqs/D:\\生信\\...\\2B3P.fa'
```
**原因**: 中文路径导致 ProteinMPNN 输出路径拼接失败
**解决**: 整个工作区搬到 `C:\Temp\mpnn_work2\` 纯英文路径

---

## 五、5 个 GFP 家族 PDB 信息

| PDB ID | 蛋白 | 长度 | 起始残基 | chromophore 位置 | 用途 |
|--------|------|------|----------|------------------|------|
| 2B3P | sfGFP | 238 | SKG (缺 M) | 64-66 (TYG) | MPNN 反折叠 |
| 2WUR | avGFP | 238 | KGE (缺 MS) | 36-38 (TYG) | MPNN 反折叠 |
| 7LG4 | amacGFP | 238 | (需检查) | (?) | 备用 |
| 2HPW | cgreGFP | 235 | (?) | (?) | 不用 (pLDDT 30) |
| 8QBJ | mBaoJin | 234 | ENM (缺 MVSK) | (?) | 不用 (chromophore 缺失) |

### 5.1 下载

```python
import urllib.request
for pdb_id in ["2B3P", "2WUR", "7LG4", "2HPW", "8QBJ"]:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    urllib.request.urlretrieve(url, f"{pdb_id}.pdb")
```

### 5.2 在 sfGFP 序列中位置

- 比赛 WT 是 sfGFP
- 长度 238 aa
- 起始: M
- chromophore: TYG (pos 65-67, 1-based)
- mBaoJin/StayGold 是完全不同的家族, Tm 92°C 来自其自身特性

---

## 六、Exclusion List 处理

### 6.1 文件位置
`D:\生信\2026Protein Design\Exclusion_List.csv`

### 6.2 规模
- 135,414 条已知 GFP 变体
- 比赛排除这些以避免"已知答案"

### 6.3 匹配规则
当前实现: **全序列完全匹配** (Python set 查询)

```python
excl = pd.read_csv("Exclusion_List.csv")
excl_seqs = set(excl["Sequence"].astype(str).str.strip())
candidate_seq in excl_seqs  # 严格匹配
```

### 6.4 未知风险
- 实际比赛可能是子串匹配 / 模糊匹配 / hash 匹配
- Round 5 验证: 用 Submission system 在线测试

### 6.5 已知要点
- mBaoJin WT 命中排除列表 (Round 3 验证)
- 需要加 1+ 突变绕开
- cgreGFP+S65T 命中排除列表 (Round 4 验证)

---

## 七、关键文献速查

### 7.1 sfGFP (Pédelacq 2006 Nat Biotechnol)

**11 个突变**:
- 折叠报告子 5: F64L, S65T, F99S, M153T, V163A
- 表面 6: S30R, Y39N, N105T, Y145F, I171V, A206V
- Tm 86.1°C (CD denaturation, BMC Res Notes 2023)
- 在 E. coli 中正确折叠

### 7.2 mBaoJin (Zhang 2024 Nat Methods)

**来源**: StayGold 单体化
**关键 8 突变** (相对 StayGold):
- S55T, H77R, E80G, Q140P, H141Q, C165Y, N171Y, T201A
- Tm 92°C
- 单体 (StayGold 是二聚体)

### 7.3 htFuncLib sf:acid.3 (Weinstein 2023 Nat Commun)

**6 突变** (相对 sfGFP):
- T65S, Q69L, S72A, T108V, Y145M, V224I
- **Tm 可达 96°C** (通过 htFuncLib 设计)
- 16,000+ 功能性变体

### 7.4 TGP (Close 2015 Proteins)

**Tm 85°C**, 90°C 半衰期 380 min
**突变** (相对 mAG):
- 稳定 7: K30I, A53S, T59P, V60A, T82A, K190E, K208R
- 表面 5: K45E, K73E, K117E, R149E, N158E

### 7.5 ProteinMPNN (Dauparas 2022 Science)

- inverse folding: 给结构, 设计序列
- 训练数据: 19700 个 native 结构
- 4 个模型权重: v_48_002, v_48_010, v_48_020, v_48_030
- v_48_020 噪声 0.2Å, 平衡选择

### 7.6 ESM-IF1 (Hsu 2022)

- 多序列比对 + 结构 → 序列
- 暂时没用到 (Round 5 候选)

### 7.7 ThermoMPNN (Dieckhaus 2024 PNAS)

- 270,000+ 训练数据
- 单点 ΔΔG 预测
- 暂时没用到 (Round 5 强烈推荐)

---

## 八、关键概念和指标

### 8.1 比赛核心公式

```
综合分 = Finit/Finit_WT × Ffinal/Finit
       = 相对亮度 × 热稳保留率
```

- Finit: CFPS 系统初始荧光 (亮度)
- Ffinal: 72°C 热处理后荧光
- Finit_WT: 同一批测量的 sfGFP WT 初始荧光

### 8.2 30% 阈值淘汰

```
若 Finit < 0.3 × Finit_WT: 该候选 0 分
```

- 即使该候选 Ffinal/Finit = 1, 0.3×1=0.3 总分 = 淘汰
- 实际影响: 极端多突变或折叠失败

### 8.3 Best Top-1 规则

6 条提交中, 取**最高分**作为队伍最终排名。意味着:
- 1 条爆款 > 6 条平均
- 策略: 集中资源到 1-2 条顶级

### 8.4 比赛奖项

- **Top-1 综合分**: 前 30% 金奖
- **最佳亮度奖**: 独立单奖项
- **最佳热稳奖**: 独立单奖项
- 独立奖项可叠加, 都需要客观分数

---

## 九、数据流和文件依赖

### 9.1 输入数据

```
D:\生信\2026Protein Design\
├── AAseqs of 5 GFP proteins_20260511.txt  # 5 个 WT 序列
├── GFP_data.xlsx                          # 14万条官方训练数据
├── Exclusion_List.csv                     # 13.5万条排除
└── work/round4/pdbs/*.pdb                 # 5 个 GFP 晶体结构
```

### 9.2 中间产物

```
work/round4/
├── candidates_round4_v3.json      # 41 手工候选
├── candidates_round4_extended.json # 53 总候选
├── esmfold_round4_v3.json         # 41 ESMFold 评估
├── esmfold_mpnn.json              # 12 sfGFP_MPNN 评估
├── esmfold_mpnn_av.json           # 8 avGFP_MPNN 评估
├── mpnn_candidates_final.json     # 12 MPNN 加载
├── mpnn_avgfp_candidates.json    # 8 MPNN 加载
└── mpnn_output_final/             # MPNN 原始 FA
```

### 9.3 最终提交

```
work/round4/
├── submission_round4_v5.csv    ⭐ 最终推荐
├── final_6_round4_v5.json      # 详细评分
├── submission_round4_v4.csv    # 备选
└── submission_round4_v3.csv    # 备选
```

---

## 十、常见问题与解答 (FAQ)

### Q1: ESMFold pLDDT 多高才"安全"?

| pLDDT | 含义 | 比赛风险 |
|-------|------|----------|
| > 80 | 极高置信, 几乎确定正确 | 极低 |
| 70-80 | 高置信, 通常正确 | 低 |
| 50-70 | 中等置信, 多数正确 | 中 (GFP 系列常在此) |
| 30-50 | 低置信, 可能未正确折叠 | 高 |
| < 30 | 极低置信, 几乎肯定错误 | 极高 |

**我们的发现**: sfGFP 系列普遍 pLDDT 40-50, 实际没问题; mBaoJin 35-40, 实际也没问题。pLDDT 是相对参考, 不是绝对阈值。

### Q2: 怎么判断一个候选的"综合比赛分"?

不能直接用 score 算, 必须用预测 Finit × Ffinal/Finit:
1. Finit = base × mut_factor × fold_prob × chromo_factor
2. Ffinal/Finit = 热稳保留率 (依赖 Tm)
3. 综合 = Finit × Ffinal/Finit

### Q3: 比赛截止日期?

未知! **强烈建议联系 root session 询问**。本轮文档基于"时间充足"假设。

### Q4: ThermoMPNN 在哪?

- GitHub: https://github.com/KuhlmanLab/ThermoMPNN
- Round 5 部署: `pip install thermompnn` 或 git clone
- 用途: 单点 ΔΔG 预测, 累加可估 multi-mutant Tm

### Q5: ESM3 API 在哪?

- https://forge.evolutionaryscale.ai
- 需要申请 token (免费 beta)
- 用途: 多模态 de novo 序列生成 (基于 ESM3)

### Q6: 比赛文档必交什么?

1. 序列文件 (CSV, 6 条, 格式严格)
2. 设计思路文档 (PDF, 流程+算法+筛选逻辑)
3. 开源仓库链接 (GitHub/GitLab/HuggingFace)
4. 根目录 README.md (环境配置+运行说明)

**状态**: Round 4 完成序列文件, **PDF 和 GitHub 待写**。

---

## 十一、附录: 关键 URL 速查

| 资源 | URL |
|------|-----|
| ESMFold (HuggingFace) | https://huggingface.co/facebook/esmfold_v1 |
| ProteinMPNN (原) | https://github.com/dauparas/ProteinMPNN |
| ProteinMPNN (国内镜像) | https://gitcode.com/gh_mirrors/pr/ProteinMPNN |
| PyTorch (国内镜像) | https://mirror.sjtu.edu.cn/pytorch-wheels/cu128 |
| RCSB PDB | https://www.rcsb.org |
| ThermoMPNN | https://github.com/KuhlmanLab/ThermoMPNN |
| ESM-IF1 | https://github.com/facebookresearch/esm |
| ESM3 Forge API | https://forge.evolutionaryscale.ai |
| htFuncLab (Weinstein 2023) | https://github.com/Fleishman-Lab |
| mBaoJin PDB (8QBJ) | https://www.rcsb.org/structure/8QBJ |

---

**完成日期**: 2026-06-22
**作者**: Trae AI Agent
