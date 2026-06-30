# Round 10 复现 Runbook

> **目标**: 从 R4 的 MPNN_T01_014 出发，迭代生成 R10 排序分 0.815 的候选

---

## 一、前置条件

- R4 最佳候选序列（MPNN_T01_014）已生成
- ProteinMPNN v_48_020 已部署到 `C:\proteinmpnn_r10`
- ESMFold (facebook/esmfold_v1) 已下载到本地
- 16GB 显存 GPU 可用

## 二、关键步骤

### 步骤 1: 生成 MPNN_T01_014 的 ESMFold 预测 PDB

```python
from transformers import AutoTokenizer, EsmForProteinFolding
import torch

tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained(
    "facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True).cuda()
model.trunk.set_chunk_size(128); model.eval()

seq = "M..."  # MPNN_T01_014 完整序列
tokens = tokenizer([seq], return_tensors="pt", add_special_tokens=False)["input_ids"].cuda()
with torch.no_grad():
    output = model(tokens)
# 保存 output.position_outputs 作为 PDB
```

### 步骤 2: 残基级 pLDDT 分析

```python
plddt_raw = output.plddt.cpu().numpy()[0]
atom_mask = output.atom37_atom_exists.cpu().numpy()[0]
plddt_scaled = plddt_raw * 100.0
masked_sums = (plddt_scaled * atom_mask).sum(axis=1)
masked_counts = atom_mask.sum(axis=1).astype(float)
plddt_per_res = masked_sums / masked_counts
# 识别 pLDDT < 70 的残基作为"低置信区域"
low_conf = [i for i, p in enumerate(plddt_per_res) if p < 70]
```

### 步骤 3: ProteinMPNN 局部重设计

```bash
cd C:\proteinmpnn_r10
python protein_mpnn_run.py \
  --pdb_path /path/to/mpnn_t01_014_pred.pdb \
  --out_folder /path/to/output \
  --num_seq_per_target 50 \
  --sampling_temp "0.1 0.2 0.3 0.5 0.7" \
  --seed 37 \
  --batch_size 25
```

### 步骤 4: ESMFold 验证（num_recycles=8）

```python
with torch.no_grad():
    output = model(tokens, num_recycles=8)
# 计算 pLDDT, chromo (残基 58-72), pTM
```

### 步骤 5: 排序分计算

```python
score = 0.40 * ptm + 0.30 * (plddt / 100) + 0.30 * (chromo / 100)
```

## 三、关键参数

| 参数 | 推荐值 | 备注 |
|:-----|:-------|:-----|
| num_recycles | 8 | 竞赛规则要求 |
| 固定残基 | T65, Y66, G67, R96, E222 | 仅 5 个核心 |
| 采样温度 | 0.1, 0.2, 0.3, 0.5, 0.7 | 多温度覆盖 |
| num_seq_per_target | 50/温度 | 每个温度 50 条 |
| 种子 | 37 | 可任意 |
| batch_size | 25 | 受显存限制 |

## 四、预期输出

- 总候选数: ~285（5 温度 × ~50-60 条）
- 生存底线通过率: ~80%
- Top 1 排序分: 0.81-0.82
- Top 6 平均排序分: 0.78-0.79

## 五、故障排查

| 问题 | 解决方案 |
|:----|:--------|
| MPNN 输出 0 条 | 检查 PDB 是否完整（238aa） |
| 全部候选 pLDDT < 60 | 减少固定残基数（→5）或降低温度 |
| 中文路径报错 | 复制到 `C:\proteinmpnn_r10` |
| ESMFold OOM | 减小 `chunk_size` 至 64 |

---

*Runbook 作者: Trae AI Agent (Claude)*
*最后更新: 2026-06-23*
