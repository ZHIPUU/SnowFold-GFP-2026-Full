# GFP 论文突变知识库 (paper_kb)

> 用途: 第二轮设计候选时,从这些论文里挑选已知能提升亮度 / 稳定性的关键突变。
> 来源: referencepaper/ 目录下 5 篇 PDF。
> 评分依据: 综合分 = Finit/Finit_WT × Ffinal/Finit(72°C 热处理后)。需要两者都高。
> **关键洞察**:训练数据集基于已部分 sfGFP 化的中间体,**纯数据驱动天花板 ~2.33× Finit/Finit_WT**。要打破第一轮 9.50,必须叠加论文知识突变。

---

## 1. Superfolder GFP (sfGFP) — Pédelacq et al. 2006

**核心信息**:
- 在 avGFP 基础上加 11 个突变得到,Tm ~78°C(WT avGFP Tm ~64°C)
- 在大肠杆菌中正确折叠,不依赖分子伴侣

**11 个突变**:
- **折叠报告子核心 5 突变**(F64L/S65T/F99S/M153T/V163A)——增强折叠和荧光
- **额外稳定性 6 突变**(S30R/Y39N/N105T/Y145F/I171V/A206V)——改善表面 loop

**关键单突变效应**(ΔG 提升):
| 突变 | ΔΔG (kcal/mol) | 备注 |
|---|---|---|
| **S30R** | **+1.25** | 显著提升热力学稳定性 |
| F64L | 中等 | 折叠关键 |
| S65T | 大 | 改造为绿色荧光,加快成熟 |
| F99S | 大 | 折叠关键 |
| M153T | 大 | 折叠关键 |
| V163A | 中等 | 折叠关键 |
| Y39N, N105T, Y145F, I171V, A206V | 小到中等 | 表面 loop 优化 |

**对设计的启示**:
1. S65T 是几乎所有"现代" GFP 的标准突变(加快成熟、提升亮度)
2. F64L/M153T/V163A 三件套是折叠关键
3. **S30R 单点就给 +1.25 kcal/mol,极有价值**
4. A206V 防止二聚化

---

## 2. TGP — Close et al. 2015 (Los Alamos)

**核心信息**:
- 由 mAG 出发,经 directed evolution 得到 eCGP123(异常稳定但聚集)
- eCGP123 加 5 个表面电荷反转突变 + C 端 GGGSGGG linker → TGP
- 85°C 几乎完全保持荧光;90°C 半衰期 ~380 min(对比 eCGP123 ~175 min, mAG 几乎瞬间失活)

**精确突变表**(从 Figure 2D alignment 直接提取,基于 mAG 编号):

### 稳定性突变(在 eCGP123 和 TGP 中均存在,7 个)
| 位置 | mAG | eCGP123/TGP | 机制 |
|---|---|---|---|
| **30** | K | **I** | β-strand 偏好,消除电荷排斥 |
| **53** | A | **S** | 扩展内部 H 键网络(Ser53 Oγ + Thr136 + Asp55) |
| **59** | T | **P** | 稳定 central helix(310-helix 二面角) |
| **60** | V | **A** | 配合 T59P,稳定 central helix |
| **82** | T | **A** | 改善 packing,减少 off-pathway folding |
| **190** | K | **E** | 表面 H 键 + 反转电荷 |
| **208** | K | **R** | 与 Asp41 形成 salt bridge |

### TGP 特有表面突变(5 个, 用于改善溶解度)
| 位置 | mAG/eCGP123 | TGP | 作用 |
|---|---|---|---|
| **45** | K | **E** | 反转电荷,消除界面接触 |
| **73** | K | **E** | 同上 |
| **117** | K | **E** | 同上 |
| **149** | R | **E** | 破坏 C-E protomer 界面(1462 Å²) |
| **158** | N | **E** | 同上 |
| **219-225** | MLPSQAK | **GGGSGGG** | C 端柔性 linker,无聚集倾向 |

**对设计的启示**:
1. **30I + 53S + 59P + 60A + 82A** 是 buried 稳定核心,适合放进所有候选
2. **190E + 208R** 是 surface H 键桥
3. **45/73/117/149/158 → E** 改善聚集(尤其在高浓度 CFPS 体系)
4. **C 端替换** 对 72°C 稳定性有贡献

---

## 3. StayGold — Hirano et al. 2022

**核心信息**:
- 来源: Cytaeis uchidae 绿色荧光蛋白
- 单体,极端光稳定性,Tm ~92°C
- 在高温、强光下保持荧光

**启示**:
- StayGold 骨架与 avGFP 同源性较低,直接突变改造可能破坏折叠
- 核心思想:可借鉴但难以直接移植到其他 GFP 骨架

---

## 4. mBaoJin — 2024 新型超稳定绿色荧光蛋白

**核心信息**:
- Tm 92°C,极强酸碱耐受
- 70°C 加热 1 小时仍保持荧光

(具体突变待补充)

---

## 5. Local fitness landscape of GFP — Sarkisyan et al. 2019

**核心信息**:
- 大规模突变扫描,展示 GFP 适应度地形
- 高适应度突变集中在特定位置
- 强 epistasis

---

## 综合设计策略

### 候选突变池(从 5 篇论文汇总)

**亮度优先**(在已有数据中可能未出现,需叠加):
- S65T (chromophore maturation, +bright)
- F64L, F99S, M153T, V163A (sfGFP 折叠报告子核心)
- I152V/T/S (chromophore 邻位,Round 1 已用)

**稳定性优先**:
- **S30R (+1.25 kcal/mol)** — 单点最猛
- **T30I** — TGP 风格(注意:S30R 和 T30I 冲突,选一个)
- **A53S** — TGP,扩展 H 键网络
- **T59P + V60A** — TGP 核心,稳定 central helix
- **T82A** — TGP,消除 off-pathway folding
- Y39N, N105T, Y145F, I171V, A206V (sfGFP 表面)
- **K190E, K208R** — TGP,表面 H 键桥

**降低聚集**:
- **K45E, K73E, K117E, R149E, N158E** (TGP 表面电荷反转)
- **C 端 GGGSGGG 替换 MLPSQAK** (TGP)

### 6 条候选的核心思路(更新版)

1. **avGFP + sfGFP 11 突变 + TGP 稳定核心**(Round 1 Seq 6 强化版)
2. **avGFP + sfGFP 折叠报告子 5 突变 + TGP 稳定核心 + 表面 E 突变**(最高亮度 + 高稳定 + 抗聚集)
3. **cgreGFP + sfGFP 折叠报告子 5 突变 + TGP 风格突变**(利用 cgreGFP 高 baseline)
4. **amacGFP + sfGFP 折叠报告子 + TGP 稳定 + S30R 替代 T30I**
5. **ppluGFP + sfGFP 折叠报告子 + TGP 稳定**
6. **avGFP + sfGFP 11 + TGP 表面突变 + C 端替换**(最强稳定性)

### 与 Round 1 差异

- Round 1 只用 sfGFP 11 个突变
- Round 2 额外加 TGP 风格突变(7 stability + 5 surface + C-term)
- 候选数:Round 1 = 6,Round 2 = 6(从数据+论文综合选)

### 风险

1. 多突变组合可能产生意外 epistasis(ESMFold 验证可降低风险)
2. 移植 sfGFP 突变到 cgreGFP/amacGFP/ppluGFP 兼容性需要 ESMFold 验证
3. TGP 风格突变基于 mAG 编号,移植到 avGFP/cgreGFP 时位置可能不严格对应(需要结构对齐)

---

## 6 条候选的具体设计(下一步生成)

设计完成后,每条都用 XGBoost 打分,选 Top-6 提交。

候选骨架(基于 ESM2-650M 嵌入 + XGBoost 排序 + 论文突变组合):
- 骨架 1: avGFP WT(基础)
- 骨架 2: cgreGFP WT(高 baseline)
- 骨架 3: amacGFP WT
- 骨架 4: ppluGFP WT
- 骨架 5: sfGFP 完整(Round 1 验证)
- 骨架 6: cgreGFP + sfGFP-style 突变(新组合)

每条骨架 + 不同论文突变组合,共生成 10-15 条候选,最后选 6。