# Round 15 踩坑记录

---

## 1. ProteinMPNN scoring 失败

### 问题
R15 计划用 ProteinMPNN 自己的 confidence 分数作为第 3 个评分维度（替代 ESM2），但 MPNN scoring 在本机失败：

```
Zero-size array in numpy.min()
```

### 根因
ProteinMPNN 的 scoring API 在 Windows + Python 3.14 + numpy 2.x 组合下报错；可能因为所有候选序列无法被 MPNN 正确 tokenize。

### 解决
用 ESM2 似然作为第 3 个评分维度（替代方案）。

### 教训
> **多模型部署耗时是隐性成本**——主备方案是必要的。MPNN scoring 与 ESM2 似然语义相近，可以互为替代。

---

## 2. ESM-IF1 模型下载超时

### 问题
原本计划部署 ESM-IF1 作为第 4 个评分维度，但 HuggingFace 下载超时（>30 min）。

### 解决
放弃 ESM-IF1，用 ESM2 似然 + ESMFold 稳定性替代。

### 教训
> 模型下载需提前规划，不能临时抱佛脚。

---

## 3. ESM2 似然与 ESMFold 不完全一致

### 问题
部分候选 ESMFold score 高但 ESM2 似然略低（反之亦然）。

| 候选 | ESMFold score | ESM2 logP |
|:-----|:-------------:|:--------:|
| R14_A_T01_013 | 0.877 | -0.3865 (略低) |
| R14_D_T02_033 | 0.869 | -0.3726 (最高) |

### 根因
ESM2 是序列似然模型，ESMFold 是结构预测模型——两者优化目标不同。

### 教训
> **共识评分比单一评分更稳健**——但共识可能掩盖"真正的"高质量候选。

---

## 4. final_score 阈值的隐患

### 问题
final_score 0.92+ 阈值在 R15 选出 6 条候选，但这阈值是"相对"的——可能在另一批候选中选出 6 条低质量序列。

### 教训
> 阈值应**绝对**而非相对（如 final_score > 0.85），而非"前 6 名"。

---

*踩坑记录作者: Trae AI Agent (Claude)*
*最后更新: 2026-06-25*
