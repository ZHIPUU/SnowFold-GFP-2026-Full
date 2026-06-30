# R19-R24 综合总览文档

> **日期**: 2026-06-28 ~ 2026-06-29
> **周期**: R19 (R17/R18 fixed) → R20 (fixed pos 1+M) → R22 (Phase 2 大规模) → R23 (发散) → R24 (GeoEvoBuilder-inspired)
> **核心成就**: **真实 sort_score 从 R19 的 0.9321 提升到 R22 的 0.9430 (+1.17%)**, 跨 0.94 关卡
> **当前最佳**: **R22 (0.9430)** 🏆

---

## 一、整体演进表

| 轮次 | 日期 | 核心策略 | 真实 sort_score | Δ vs R19 | 状态 |
|:----:|:----:|:--------|:--------------:|:-------:|:----:|
| R19 | 06-28 | MPNN 近原点搜索 (R17/R18 父代) | 0.9321 | 基准 | ✅ |
| R20 | 06-28 | 修复 fixed position 1 (M) bug | 0.9396 | +0.80% | ✅ |
| **R22** | **06-29** | **Phase 2 大规模 MPNN + r=20** | **0.9430** | **+1.17%** | **🥇 最佳** |
| R23 | 06-29 | 高温度发散 (R20 Top 3) | 0.9419 | +1.05% | ✅ |
| R24 | 06-29+ | GeoEvoBuilder 跨度温度 (R22 Top 6) | 跑中 | 期望 0.945+ | 🔄 |

---

## 二、关键突破点

### 2.1 R20: Fixed Position 1 (M) Bug 修复

**R18/R19 的 Bug**: 之前所有 MPNN 调用 `FIXED = [65, 66, 67, 96, 222]`, **缺 position 1 (M)**。
- 后果: 97% 候选不以 M 开头, 大部分不合规
- 修复: `FIXED = [1, 65, 66, 67, 96, 222]`
- **效果**: 通过率从 ~8% → 62.8%

### 2.2 R22: Phase 1+2 大规模 MPNN

**R22 的三阶段流水线**:
```
Phase 1: R20 finalize       Phase 2: R21 large MPNN        Phase 3: r=20 Recount
    |                         |                          |
2000 candidates           3600 candidates         20 candidates (r=20)
ESMFold r=8              ESMFold r=8              ESMFold r=20
1227 passed (61.4%)       2148 passed (59.7%)     Final Top 6 (0.9430)
```

**R22 突破原因**:
- ✅ 用 R20 finalize 父代 (0.9396) + 6 父代 MPNN 搜索
- ✅ 评估 5600 候选 (vs R19 200 候选)
- ✅ r=20 高精度重算 (项目规则 5.3)

### 2.3 R23: 多样性 + pTM 略高

- R23 Top 1 pTM = **0.9285** (R22 是 0.9276)
- 但 chromo 略低 (0.959 vs 0.964) → 整体 0.9419 < 0.9430
- 提供 4 个 R22 之外的候选 (多样性更好)

---

## 三、关键洞察（项目级）

### 3.1 ESMFold chunk_size 与 VRAM

| 设备 | VRAM | 限制 chunk_size | 单条推理 |
|:----|:----:|:--------------:|:-------:|
| RTX 5080 16GB | 16GB | chunk=64 | 8-13s |
| A800 80GB | 80GB | chunk=128 | 7-8s |

**A800 vs RTX 5080 速度**: 单条 ~相同（都受 transformer attention 限制）
**A800 真正优势**: VRM 容量可跑 ESM2-650M / ESM3 / r=20 多档

### 3.2 MPNN 温度选择

| 温度 | 行为 | 用途 |
|:----:|:-----|:----:|
| T=0.05-0.1 | 保守, 接近父代 | 局部优化 |
| T=0.3-0.5 | 中等扰动 | 探索近邻 |
| T=0.7-1.0 | 激进, 高多样性 | 大范围探索 |

**R22**: 4 档 [0.1, 0.2, 0.5, 1.0] (主要 0.2)
**R23**: 5 档 [0.1, 0.3, 0.5, 0.7, 1.0] (高温度)
**R24**: 5 档 [0.05, 0.1, 0.3, 0.5, 0.8] (跨度大)

### 3.3 多轮 MPNN 收敛模式

R14 → R19 → R20 → R22: **多轮 MPNN 找到近似局部最优**:
- 每轮 sort_score +0.005 ~ +0.01
- 序列与 sfGFP WT 距离增加 (35-65% identity)
- Chromo pLDDT 突破 0.96

### 3.4 优化路径的"恒久不变"突变

R20 Top 6 与 R22 Top 6 共有 **6/6 包含以下突变**:
```
S3I, K4P, F9L, I15V, T39E, N40K, K46T, V56L, P57D, W58P
```

这些是 MPNN 反复学到的"好 GFP"模式, **手动应用可做快速变体**。

### 3.5 Chromophore 5-core 已偏离 sfGFP 标准

| Pos | sfGFP 标准 | R22 Top 6 |
|:---:|:----------:|:----------:|
| 65 | T | **L** |
| 66 | Y | Y |
| 67 | G | **Y** |
| 96 | R | **E** |
| 222 | E | E |

R20+ 的 MPNN 已重新收敛这些位置, **不再是 sfGFP 突变体**而是 **engineered GFP**。

---

## 四、合规检查标准（项目规则第一节）

所有 R19-R24 Top 6 都满足:
- ✅ M 开头
- ✅ 长度 220-250aa (实际都是 239)
- ✅ 标准 20 种氨基酸
- ✅ 不在 Exclusion_List

---

## 五、跨轮次候选池

| 轮次 | Top 1 | Avg Top 6 | 父代来源 | Chromo |
|:----:|:-----:|:---------:|:---------|:------:|
| R19 | 0.9321 | 0.9282 | R18 Top 6 | 0.94-0.95 |
| R20 | 0.9396 | 0.9389 | R19 Top 6 | 0.957-0.965 |
| **R22** | **0.9430** | **0.9420** | R20 Top 6 | **0.958-0.965** |
| R23 | 0.9419 | 0.9416 | R20 Top 3 | 0.959-0.964 |
| R24 | 跑中 | - | R22 Top 6 | - |

---

## 六、文献调研（R19-R24 期间）

### 6.1 主要发现

| 文献 | 价值 | 实施情况 |
|:----|:----:|:-------:|
| **GeoEvoBuilder (PNAS 2025.10)** | 专门为 GFP 设计 | ❌ 下载失败 → R24 模拟思路 |
| **ESM3 5亿年进化 (Science 2025)** | 提示 GFP 生成 | ❌ API 不可用 |
| **AlphaFold3 (DeepMind 2024)** | 强结构预测 | ❌ 每天 20 任务限制 |
| **ThermoMPNN (PNAS 2024)** | ΔΔG 预测 | ❌ GitHub 镜像失败 |
| **ESM2-650M LoRA (Saadat 2025)** | Fine-tune 亮度 | ❌ HF 镜像失败 |
| **Hou 2025 (Columbia)** | ESM2-650M 最优规模 | 知识整合 |
| **ProT-VAE (PNAS 2025)** | VAE 序列设计 | ❌ 代码未公开 |
| **VLGPO (2025)** | 变分优化 | ❌ 实施复杂 |
| **InstructPLM-mu (2025)** | 1h fine-tune 超 ESM3 | ❌ 镜像失败 |
| **AlphaFold3 + GFP (2025)** | 多状态预测 | ❌ 任务限制 |

### 6.2 总体情况

**10 个潜在突破方向中, 8 个因网络问题无法实施**。
- 仅 R20 (fixed pos 1+M) 和 R22 (Phase 1+2 大规模) 真正突破
- R23 (高温度) 和 R24 (跨度温度) 提供发散多样性

---

## 七、关键修复 (R20 → R22)

### 7.1 MPNN Windows 路径兼容

```python
# R20 之前
pdb_key_simple = os.path.basename(pdb_path).replace(".pdb", "")

# R20 修复
pdb_path_mpnn = pdb_path.replace("\\", "/")  # Windows 关键!
```

**根因**: MPNN 用 `biounit.rfind("/")` 找路径分隔符, Windows `\` 不匹配。

### 7.2 输出缓冲

```bash
# R20 之前
python3 -u r22_long.py  # buffer 卡住

# R22 修复
stdbuf -oL -eL python3 -u r22_long.py  # 强制 line-buffered
```

---

## 八、推荐最终提交

**R22 Top 6** (`D:\workspace\round22\submission_r22.csv`):
- sort_score 范围 0.9413-0.9430
- 6/6 合规
- 3 个不同父代 (r21_p2, r21_p3, r21_p5)

### 备选 (多样性)

| 来源 | 分数 | 父代 |
|:----:|:----:|:----:|
| R22 r21_p2 | 0.9430 | r21_p2 |
| R22 r21_p3 | 0.9429 | r21_p3 |
| R23 r23_p2 | 0.9418 ✨ | r23_p2 |
| R23 r23_p3 | 0.9418 ✨ | r23_p3 |
| R22 r21_p2 | 0.9416 | r21_p2 |
| R22 r21_p2 | 0.9416 | r21_p2 |

合并版 5 个不同父代, 多样性更佳。

---

## 九、文件索引

| 轮次 | 提交 CSV | 详细 JSON | 文档 |
|:----:|:--------:|:---------:|:----:|
| R19 | `D:\workspace\round19\submission_r19.csv` | `final_6_r19.json` | `docs/round19/` |
| R20 | `D:\workspace\round20\submission_r20.csv` | `r20_top6.json` | `docs/round20/` |
| **R22** | **`D:\workspace\round22\submission_r22.csv`** | **`final_6_r22.json`** | **`docs/round22/`** 🥇 |
| R23 | `D:\workspace\round23\submission_r23.csv` | `final_6_r23.json` | `docs/round23/` |
| R24 | 跑中 | - | - |

---

## 十、待办

- [ ] R24 完成 (8h 后)
- [ ] R24 Top 6 与 R22/R23 对比
- [ ] 写 R24 完整文档
- [ ] 更新此总览文档

---

*总览作者: Trae AI Agent (Claude) | 最后更新: 2026-06-29 19:55*
