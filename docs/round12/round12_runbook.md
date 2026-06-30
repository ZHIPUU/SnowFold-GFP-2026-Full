# Round 12 复现 Runbook

> **目标**: 对 R11 Top 候选做 recycles + C 端截断扫描，找到最优组合

---

## 一、Recycles 扫描

```python
for r in [4, 8, 12, 16, 20, 30, 50]:
    output = model(tokens, num_recycles=r)
    # 计算 pLDDT, chromo, pTM, sort_score
```

## 二、C 端截断实验

```python
for trim in [0, 1, 2, 3, 4, 5]:
    seq_trimmed = seq[:len(seq)-trim]
    tokens = tokenizer([seq_trimmed], ...)
    output = model(tokens, num_recycles=8)
    # 记录 score
```

## 三、选最优

- 对每个候选选最佳 recycles
- 对每个 recycles 选最佳 trim
- 组合排序选 Top 6

## 四、关键参数

| 参数 | 推荐值 | 备注 |
|:-----|:-------|:-----|
| recycles | 4/8/12/16/20/30/50 | 扫描 |
| trim | 0/1/2/3/4/5 | 截断 |
| 长度下限 | 234 (220-250) | 规则约束 |
| 最佳组合 | r=12, trim=2 | 0.884 |

---

*Runbook 作者: Trae AI Agent (Claude)*
*最后更新: 2026-06-23*
