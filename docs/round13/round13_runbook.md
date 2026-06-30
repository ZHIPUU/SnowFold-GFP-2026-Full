# Round 13 复现 Runbook

> **目标**: 通过"截断扫描 + MPNN 新种子"双路径策略保证分数 + 提升多样性

---

## 一、路径 A: 截断扫描

```python
# 对 R11 Top 10 候选做 trim
for candidate in r11_top10:
    for trim in [0, 1, 2, 3, 4]:
        seq_trimmed = candidate["seq"][:len(candidate["seq"])-trim]
        # ESMFold 评估
        score = evaluate(seq_trimmed)
```

## 二、路径 B: MPNN 新种子

```bash
python protein_mpnn_run.py \
  --pdb_path r10_t02_078_pred.pdb \
  --out_folder mpnn_output_seed137/ \
  --num_seq_per_target 50 \
  --sampling_temp "0.1 0.2 0.3 0.5" \
  --seed 137 \
  --fixed_positions_jsonl fixed_5core.jsonl
```

## 三、混合策略

- 4 条来自路径 A（最高分）
- 2 条来自路径 B（独立验证）

## 四、关键参数

| 参数 | 推荐值 | 备注 |
|:-----|:-------|:-----|
| trim 范围 | 0-4 | 局部最优在 trim2 |
| seed | 137 | 与 R10 seed=37 区分 |
| 温度 | 0.1-0.5 | 多温度 |
| 多样性比例 | 4:2 (A:B) | 30%+ 不同源 |

---

*Runbook 作者: Trae AI Agent (Claude)*
*最后更新: 2026-06-23*
