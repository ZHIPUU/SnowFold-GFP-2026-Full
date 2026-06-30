# Round 11 复现 Runbook

> **目标**: 从 R10 最佳候选出发，迭代生成 R11 排序分 0.881 的候选

---

## 一、关键步骤

### 步骤 1: 准备 R10 最佳候选

```python
# R10r3_T02_078
seq = "MPIPGDELLSGVVPVKVNLDGDVNGNKFKIKGEGTGDATKGELKLTFKVTEGELPLDWVLIVDILTYGLRIFWKLPEDNPLRDFYKACLPEGYKIERTLKFKDEGTLTVTSDVRFEGDTLVSDIELKGTDFKEGGLLLGKTVATLTYSGKVEVSPDEEKHGVKLTYTLTLTFTDGSTATADVEELISPIGKGPETLPEPRTLYLRRVYSDDPTAKGDTVVVEQKQEEGVPLPPPELLEA"
```

### 步骤 2: ESMFold 预测 R10 结构

```python
output = model(tokens, num_recycles=8)
plddt_per_res = ...  # 残基级 pLDDT
# 固定 pLDDT >= 70 的 189 个残基 + 5 核心
fixed_positions = [i for i, p in enumerate(plddt_per_res) if p >= 70]
```

### 步骤 3: 生成 fixed_positions JSONL

```json
{"label": "fixed", "positions": [64, 65, 66, 95, 221, ...]}
```

**注意**: 索引是 0-based，5 个核心残基（T65/Y66/G67/R96/E222）对应索引 [64, 65, 66, 95, 221]。

### 步骤 4: ProteinMPNN 多温度采样

```bash
python protein_mpnn_run.py \
  --pdb_path r10_t02_078_pred.pdb \
  --out_folder mpnn_output/ \
  --num_seq_per_target 50 \
  --sampling_temp "0.1 0.2 0.3 0.5" \
  --seed 37 \
  --fixed_positions_jsonl fixed_positions.jsonl
```

### 步骤 5: ESMFold 验证 num_recycles=8

对 147 候选全部运行 ESMFold r=8，记录 pLDDT, chromo, pTM。

---

## 二、关键参数

| 参数 | 推荐值 | 备注 |
|:-----|:-------|:-----|
| num_recycles | 8 | 竞赛规则 |
| 固定策略 | pLDDT≥70 + 5 核心 | 189 个固定 |
| 采样温度 | 0.1, 0.2, 0.3, 0.5 | T=0.3 最佳 |
| num_seq_per_target | 50/温度 | 共 ~150-200 |
| 种子 | 37 | 可任意 |

## 三、预期输出

- 总候选: ~150
- 生存底线通过率: 95%+
- Top 1 排序分: 0.86-0.88
- Top 6 平均: 0.85-0.86

## 四、故障排查

| 问题 | 解决方案 |
|:----|:--------|
| 候选 pLDDT < 60 | 减少固定残基数（→ 5 核心） |
| 序列截短 | 检查 PDB 是否完整 238aa |
| 排序分 < 0.80 | 调整温度范围至 0.2-0.3 |

---

*Runbook 作者: Trae AI Agent (Claude)*
*最后更新: 2026-06-23*
