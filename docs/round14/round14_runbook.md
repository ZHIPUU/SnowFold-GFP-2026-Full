# Round 14 复现 Runbook

> **目标**: 4 路并行探索（截断/Recycles/Chromophore/MPNN 新种子）

---

## 一、路径 A: 截断 + MPNN 重新设计

```python
# 用 R11p1_T03_135_trim2 序列作为新父代
parent_seq = "MMIPGDELLS..."  # R11p1_T03_135_trim2
parent_pdb = esmfold_predict(parent_seq)  # ESMFold 预测

# MPNN 重新设计
mpnn_run(parent_pdb, fixed_5core, seed=37, temps=[0.1,0.2,0.3])
```

## 二、路径 B: Recycles 测试

```python
for r in [4, 8, 12, 16, 20, 30]:
    output = model(tokens, num_recycles=r)
    # 记录 score
```

## 三、路径 C: Chromophore 微调

```python
# 固定核心 + 重设计生色团周围 30-40 个残基
chromo_neighbors = list(range(60, 80)) + list(range(90, 100)) + list(range(220, 230))
fixed = [64, 65, 66, 95, 221]  # 5 核心
designable = [i for i in chromo_neighbors if i not in fixed]
mpnn_run(parent_pdb, fixed=fixed, designable=designable, seed=37)
```

## 四、路径 D: MPNN 新种子

```python
# 用 R10 最佳骨架 + 新 seed
parent_pdb = "r10_t02_078_pred.pdb"
for seed in [137, 237, 337]:
    mpnn_run(parent_pdb, fixed_5core, seed=seed, temps=[0.1,0.2,0.3])
```

## 五、关键参数

| 路径 | seed | 温度 | 父代 | 候选数 | 最佳 score |
|:----:|:----:|:----:|:----:|:------:|:---------:|
| A | 37 | 0.1-0.3 | R11 trim2 | ~100 | 0.892 |
| B | — | — | R11 trim2 | ~30 | 0.886 |
| C | 37 | 0.1-0.3 | R11 trim2 | ~50 | 0.870 |
| D | 137/237/337 | 0.1-0.3 | R10 | ~100 | 0.888 |

---

*Runbook 作者: Trae AI Agent (Claude)*
*最后更新: 2026-06-24*
