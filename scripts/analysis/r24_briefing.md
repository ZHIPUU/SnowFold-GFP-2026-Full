# R24 Briefing - GeoEvoBuilder-Inspired Exploration

> **日期**: 2026-06-29
> **动机**: GeoEvoBuilder (PKU 2025 PNAS) 思路 - zero-shot 同时优化热稳定性+活性
> **代码问题**: GeoEvoBuilder GitHub (https://github.com/PKUliujl/GeoEvoBuilder) 下载失败 (网络限制)
> **方案**: 复现 GeoEvoBuilder 核心思想,用现有工具 (ESMFold + MPNN)

---

## 一、GeoEvoBuilder 核心思想 (Liu et al. PNAS 2025.10)

[论文链接](https://www.pnas.org/doi/10.1073/pnas.2504117122) | [代码](https://github.com/PKUliujl/GeoEvoBuilder)

**核心思路**:
1. **结构 + 进化双约束** — adaptively integrates structural and evolutionary constraints
2. **zero-shot** — 不需要迭代实验循环
3. **同时优化活性 + 热稳定性** — 在 GFP 任务上**已经实验验证成功**
4. **支持 30% 残基变化** — 探索巨大序列空间

**对项目最直接的价值**:
- 我们的项目评分 = `Finitial/FinitialWT × Ffinal/Finitial`
  - 需要高亮度
  - 需要高热稳定性
- **与 GeoEvoBuilder 的设计目标 100% 重合!**

---

## 二、网络问题与备选方案

### 失败的下载尝试

| 镜像 | 结果 |
|:----|:----|
| `git clone https://gh-proxy.polaris-lab.com/...` | ❌ "could not determine hash algorithm" |
| `curl https://github.com/...main.zip` | ❌ "Connection timeout after 30000 ms" |
| `wget https://gh-proxy.imxbt.com/...` | ❌ 0 字节 |
| 之前: `ghproxy.com`, `hf-mirror.com`, `mirror.ghproxy.com` | ❌ 全部失败 |

### 关键决定: 基于现有工具模拟 GeoEvoBuilder 思路

我们已经有 **ESMFold + ProteinMPNN + R22 父代 (0.9430)** — 本质上 GeoEvoBuilder 思路的子集。

**R24 复现**:
- 父代: R22 Top 6 (sort_score 0.9430) — 项目最高分父代
- 温度: [0.05, 0.1, 0.3, 0.5, 0.8] (5 档, 跨度比 R22/R23 大)
- 固定: [1, 65, 66, 67, 96, 222] (含 M + 5 chromophore)
- 候选: 5 温度 × 200 = 1000/父代 × 6 父代 = **6000 候选**
- Phase 3: Top 20 用 r=20 重算
- **预计 ~7-8 小时完成**

---

## 三、为什么 R24 有希望突破 0.95

**理由**:
1. **R22 父代质量** (0.9430) > R20 父代 (0.9396) — 起点更高
2. **温度跨度大** (0.05-0.8) — 探索 R22/R23 未覆盖的温度空间
3. **多样化探索** — 6 父代 × 5 温度 = 30 个 (父代, 温度) 组合

**潜在提升**:
- 0.9430 → 0.9450+ (Δ +0.2%) 概率中等
- 0.9430 → 0.9500+ 概率较低 (需要新发现)

---

## 四、其他调研发现

### AlphaFold3 (DeepMind 2024)
- **免费** (非商业研究)
- **限制**: 每天 20 个任务/IP — **不适合批量 GFP 设计** (我们一个 round 就 ~6000 候选)
- 不考虑

### VLGPO (Bogensperger 2025)
- 变分生成式蛋白优化
- 状态空间 VAE + flow matching
- 与本项目 MPNN/ESMFold 路线不同 — 实现复杂
- 不考虑

### ProT-VAE (Sevgen 2025 PNAS)
- Transformer + VAE
- 已实验验证 PAH 酶活 2.5×, γ-carbonic anhydrase +61°C Tm
- 代码未公开
- 不考虑

### AlphaFold3 重计算
- 比 ESMFold 更准确 (DeepMind 自家 benchmark)
- 但每天 20 个任务限制
- 不考虑

---

## 五、预期时间表

| 时间 | 任务 |
|:----|:----|
| ~18:00 | R23 跑完 (R20 Top 3 父代, 高温度) |
| ~22:00 | R24 跑完 (R22 Top 6 父代, 跨度温度) |
| ~22:00 | R22 仍为最佳 (sort_score 0.9430) |
| 提交 | 用 R22 Top 6 (`D:\workspace\round22\submission_r22.csv`) |

---

## 六、文件

- `D:\workspace\r24_server.py` - R24 脚本 (服务器版)
- 远程: `/root/autodl-tmp/r24_server.py`
- 远程: `/tmp/r24.log` (实时进度)
- 任务 ID: `15e464541126`

---

*分析: Trae AI Agent (Claude) | 时间: 2026-06-29 16:00*
