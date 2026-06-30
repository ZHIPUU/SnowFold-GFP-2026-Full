# Round 4 下一步方向与战略指南

> **目的**: 给 Round 5 及后续 AI 接手者明确的优先级和决策依据
> **适用**: 改进 Round 4 提交 / 准备下一轮 / 战略决策

---

## 零、优先级总览

### 0.1 ROI 矩阵

| 任务 | 难度 | 时间 | 预期收益 | 优先级 |
|------|------|------|----------|--------|
| **设计 PDF** | 中 | 1-2h | 比赛必交 | ⭐⭐⭐⭐⭐ |
| **GitHub 仓库** | 中 | 30min | 比赛必交 | ⭐⭐⭐⭐⭐ |
| **联系 root 询问 CFPS 协议** | 低 | 5min | 高 (信息) | ⭐⭐⭐⭐⭐ |
| **部署 ThermoMPNN** | 中 | 1-2h | 中 (更准 Tm) | ⭐⭐⭐⭐ |
| **跑 MPNN T=0.5, T=0.7** | 中 | 30-60min | 中 (更多高分) | ⭐⭐⭐⭐ |
| **PROSS 在线服务** | 低 | 30min | 中 (5稳定化设计) | ⭐⭐⭐ |
| **ESM-IF1 inverse folding** | 中 | 2-3h | 中 | ⭐⭐⭐ |
| **ColabFold 验证 mBaoJin** | 中 | 1h | 低 (已是 Top-6) | ⭐⭐ |
| **ESM3 Forge API** | 高 | 1-2h | 中 (de novo) | ⭐⭐⭐ |
| **DCI_asym GNN** | 高 | 3-5 天 | 低 (Round 2 已证) | ⭐ |
| **多任务学习** | 高 | 2-3 天 | 中 | ⭐⭐ |

### 0.2 时间预算策略

| 剩余时间 | 推荐优先级 |
|---------|-----------|
| **< 6 小时** | PDF + GitHub 必交, **不改 v5** |
| **6-24 小时** | + 部署 ThermoMPNN, 重估 Tm, 可能调整 v5 |
| **1-3 天** | + 跑更多 MPNN, PROSS, ESM-IF1 |
| **3-7 天** | 全面重做 Round 5 (新策略) |
| **> 7 天** | 尝试所有未做工具, 选最优 |

---

## 一、立即执行 (必做, < 30 分钟)

### 1.1 联系 root 询问比赛细节

**为什么**: 比赛规则文档不完整, 关键参数未明。

**应问问题** (5 个最关键):

1. **比赛截止日期?**
   - 影响: 决定剩余时间预算
   - 已知: 不明 (假设充足)

2. **CFPS 系统的具体条件?**
   - 温度: 室温 / 30°C / 37°C?
   - 反应时间: 1h / 2h / overnight?
   - 缓冲: 已知 lysate 还是特定配方?
   - 影响: Finit 测量时间点, 影响 chromophore 成熟

3. **72°C 热处理时长?**
   - 30 min / 1 h / 2 h?
   - 影响: Ffinal/Finit 估值
   - **30 min** 假设下, sfGFP 保留 ~85% (Tm-T = 14°C)
   - **2 h** 假设下, sfGFP 保留 ~50% (接近 Tm)
   - 这个数据 **改变** 几乎所有候选的预测综合分

4. **比赛打分公式细节?**
   - 30% 阈值淘汰是硬阈值还是软?
   - Finit 测量是在 chromophore 完全成熟后, 还是固定时间?
   - 影响: 多突变候选的 Finit 估值

5. **已提交的 Round 4 v5 是否被接受?**
   - 还能改吗?
   - 还是需要 Round 5?
   - 影响: 决定工作方向

**问询脚本**:
```
"Hi root, 我有几个关于比赛规则的问题:
1. 截止日期是什么时候?
2. Ffinal 加热 72°C 持续多久?
3. CFPS 反应在什么温度、时间下进行?
4. 30% 阈值淘汰是 hard 阈值还是 soft 阈值?
5. 我们现在提交的 6 条 (work/round4/submission_round4_v5.csv) 是否最终, 还能改吗?
谢谢!"
```

### 1.2 写设计思路 PDF (1-2 小时)

**比赛要求** (来自 `work/phase5/design_doc.pdf` 模板):

PDF 必须包含:

1. **任务理解与目标拆解** (10%)
   - 综合分 = Finit/Finit_WT × Ffinal/Finit
   - 极低亮度阈值 (0.3×WT) 淘汰规则
   - 奖项: Top-1 + 最佳亮度 + 最佳热稳

2. **算法管线** (25%)
   - 输入数据 (14万条官方数据 + 5 个 WT + 5 个 PDB)
   - 候选生成 (文献驱动 + 手工 + ProteinMPNN)
   - 评估 (ESMFold pLDDT/pTM + 文献 Tm 估值)
   - 评分函数 (5 维度加权)
   - Top-6 选择 (多样性约束)

3. **如何平衡亮度和稳定性目标** (15%)
   - 文献先验: sfGFP (Tm 86°C, baseline 1.0) + mBaoJin (Tm 92°C, baseline 1.2) + htFuncLib (Tm 96°C)
   - 多机制叠加: htFuncLib + S30R + I152S + Q69L + S72A
   - 经验: 突变数 ≤ 13 (Arcadia 2025)

4. **筛选 6 候选的原因** (20%)
   - 详细解释每条 (见 docs/round4/round4_design_rationale.md)
   - 表格: 序列、pLDDT、Tm、突变数、机制

5. **LLM Agent 逻辑树 + 关键日志** (30%) ⭐
   - 这是比赛**最看重**的部分
   - 必须详细展示:
     - 整个思考流程
     - 关键决策点 (v1→v2→v3→v4→v5)
     - 失败教训 (Round 2 ML 失效, 中文路径坑)
     - 工具组合策略
   - 引用代码片段 + 关键输出

**PDF 生成建议**:
- 用 reportlab (Python)
- 字体: 中文用 SimSun/微软雅黑
- 页数: 5-10 页 (不要太少或太多)
- 文件名: `YourTeamName_Design_Doc_Round4.pdf`

**模板结构** (从 Round 1 借鉴):

```python
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4

doc = SimpleDocTemplate("YourTeamName_Design_Doc_Round4.pdf", pagesize=A4)
styles = getSampleStyleSheet()
story = []

# 封面
story.append(Paragraph("2026 Protein Design — Round 4 Design Document", styles["Title"]))
story.append(Spacer(1, 0.5*inch))
story.append(Paragraph("Team: YourTeamName", styles["Normal"]))
story.append(Paragraph("Date: 2026-06-22", styles["Normal"]))
story.append(Spacer(1, 0.3*inch))

# 1. 任务理解
story.append(Paragraph("1. 任务理解与目标拆解", styles["Heading1"]))
# ... 5-10 段

# 2. 方法论
story.append(Paragraph("2. 算法管线", styles["Heading1"]))
# 数据 / 加性 / ESMFold / MPNN / 评分

# 3. 平衡
story.append(Paragraph("3. 亮度和稳定性的平衡", styles["Heading1"]))
# ...

# 4. 6 候选详解 (含表格)
story.append(Paragraph("4. 最终 6 候选详解", styles["Heading1"]))
data = [
    ["Seq", "Name", "Scaffold", "Mut", "pLDDT", "Tm", "Notes"],
    ["1", "MPNN_T01_014", "sfGFP_MPNN", "57", "68.3", "88", "ProteinMPNN de novo 王者"],
    # ...
]
story.append(Table(data))
# 每条详细解释

# 5. LLM Agent 逻辑树 (最重要!)
story.append(Paragraph("5. LLM Agent 逻辑树与执行日志", styles["Heading1"]))
# v1 → v2 → v3 → v4 → v5 演进过程
# 关键决策截图
# 失败教训
# 工具组合

doc.build(story)
```

### 1.3 建立 GitHub 仓库 (30 分钟)

**必含文件**:
- `README.md` - 项目说明 + 环境 + 运行
- `docs/round4/*.md` - 6 份文档
- `work/round4/submission_round4_v5.csv` - 最终提交
- 关键脚本 (01/06/12/19/21)
- `requirements.txt`
- `LICENSE` (MIT)

**README.md 模板** (已详细写在 runbook.md § 7.3)

**关键**:
- 公开仓库 (比赛要求"可以复现")
- 根目录 README.md (比赛明确要求)
- 包含完整运行说明

---

## 二、短期改进 (1-3 天)

### 2.1 部署 ThermoMPNN (1-2 小时)

**目的**: 用真实 ΔΔG 预测替代 v5 的 pLDDT 间接估值 Tm

**步骤**:
```bash
# 1. 克隆
git clone https://github.com/Kuhlman-Lab/ThermoMPNN.git
cd ThermoMPNN
pip install -e .

# 2. 下载模型权重
# (通常 200-500MB, 国内镜像可能需要备用)

# 3. 写预测脚本
python -c "
from thermompnn import ThermoMPNN
model = ThermoMPNN()
# 对 v5 Top-6 跑全 238 个位点的单点 ΔΔG
# 累加得到 multi-mutant ΔΔG
# 用 ΔΔG → ΔTm 公式
"
```

**预期输出**:
- 每条候选的精确 ΔΔG 和 ΔTm
- 可能改变 v5 的最终 Top-6 排序

**风险**:
- ThermoMPNN 模型大 (500MB+), 下载可能慢
- 国内可能没有镜像
- Round 5 时间预算紧张时, 跳过

### 2.2 跑更多 MPNN 温度 (30-60 分钟)

**目的**: 用 T=0.5 (激进入化) 和 T=0.7 (高多样性) 探索更激进的 de novo 设计

**实施**:
```python
# 修改 09e_proteinmpnn_v5.py
configs = [("0.1", 30, "T01"), ("0.3", 30, "T03"), ("0.5", 30, "T05"), ("0.7", 20, "T07")]
# T=0.5 应该产生更多"突变密集"但可能 pLDDT 较低的候选
# T=0.7 产生高度多样但折叠不可靠
```

**预期**:
- 50-80 条新 MPNN 候选
- 评估后可能找到 pLDDT 65+ 的新候选
- 替换 v5 中的次优候选

### 2.3 PROSS 在线服务 (30 分钟)

**目的**: Rosetta 物理能量函数 + 进化保守性约束 = 多个稳定化设计

**步骤**:
1. 访问 https://pross.weizmann.ac.il/
2. 上传 sfGFP 2B3P PDB
3. 选择 5-10 个输出
4. 下载 FASTA
5. 评估 + 加入候选池

**预期**:
- 5 个稳定化变体 (每个 8-30 突变)
- Tm 估计提升 5-15°C
- 可能提供 mBaoJin/TGP 类的高 Tm 候选

**优势**:
- 0 安装成本
- 物理级别稳定化 (Rosetta)
- 与文献方法互补

### 2.4 ESM-IF1 inverse folding 打分 (2-3 小时)

**目的**: 给每条候选打"逆向折叠概率"分数 (ESM-IF1 给定 backbone, 设计序列)

**步骤**:
```python
# ESM-IF1 安装
git clone https://github.com/facebookresearch/esm.git
cd esm
pip install -e .

# 加载模型
from esm.inverse_folding import MultimerStructureDecoder
decoder = MultimerStructureDecoder("esm_if1_gvp4_t12_100M_t2_tev_saufif_20220404")

# 对每条候选: 给出 backbone, 问"该序列的逆折叠似然"
score = decoder.score(seq, structure_coords)
```

**预期**:
- 每条候选的 inverse folding log-likelihood
- 这是第 2 维度的"折叠正确性"指标
- 与 ESMFold pLDDT 互补 (pLDDT 看序列预测的结构, IF1 看结构对应的序列)

---

## 三、中期改进 (3-7 天, Round 5 重做)

### 3.1 重新设计 Round 5 策略

**保留 v5 的核心**:
- MPNN de novo 路线
- 多骨架多样性
- ESMFold pLDDT 验证
- 不用 ML 打分

**改进点**:

1. **跑 200+ MPNN 候选** (vs 现在的 20)
   - 选 Top 6 而非 Top 6 from 20
   - 期望找到 2-3 条 pLDDT 65+ 的"新王"

2. **加入 PROSS / ProteinMPNN-Tm 联合**
   - PROSS 5 个稳定化变体
   - MPNN 重设计 100+ 候选
   - 用 ESMFold 排序

3. **加入 amacGFP MPNN** (Round 4 跳过)
   - 7LG4 PDB chain A 长度 238
   - 用 ProteinMPNN 反折叠
   - 期望 pLDDT 40-55

4. **使用 ThermoMPNN 真实 ΔΔG**
   - 替代 v5 间接 pLDDT-Tm 估值
   - 给出每条候选的精确 Tm 范围

5. **v6 重排名**
   - 重新选 Top-6
   - 综合 MPNN + PROSS + 手工

### 3.2 评估 4 个 PNNN de novo 骨架

| 骨架 | 已有 MPNN | 期望 pLDDT | 状态 |
|------|----------|----------|------|
| sfGFP 2B3P | ✅ 12 条 | 60-70 | 完成 |
| avGFP 2WUR | ✅ 8 条 | 55-65 | 完成 |
| amacGFP 7LG4 | ❌ | 40-55 | Round 5 候选 |
| cgreGFP 2HPW | ❌ | 30-40 | 跳过 (pLDDT 太低) |
| mBaoJin 8QBJ | ❌ | N/A | chromophore 缺失 |

**目标**: Round 5 跑 amacGFP MPNN, 探索第 4 个 MPNN 骨架

### 3.3 引入 ESM3 API (2-3 小时)

**目的**: ESM3 是最新最强多模态 de novo 蛋白设计

**步骤**:
1. 申请 Forge API token (https://forge.evolutionaryscale.ai)
2. 写调用脚本

```python
from esm.sdk import client
from esm.sdk.api import ESMProtein, GenerationConfig

model = client(model="esm3-medium-2024-03", url="https://forge.evolutionaryscale.ai", token=token)

# 固定 chromophore TYG
prompt_seq = "_" * 237
prompt_seq[64] = "T"  # 65 - 1 = 64 (0-based)
prompt_seq[65] = "Y"
prompt_seq[66] = "G"

protein = ESMProtein(sequence=prompt_seq)
# 关键结构位
protein.structure_predictions_enabled = True
```

3. 调用 ESM3 生成
4. 评估 + 加入候选池

**风险**:
- API 限制 (可能收费)
- 中文网络可能连不上
- Round 5 时间紧张时跳过

---

## 四、长期探索 (1-2 周, Round 6+)

### 4.1 DCI_asym GNN (3-5 天)

**目的**: 用动力学驱动的 epistasis 模型 (Round 2 已证 ML 失效, 但 DCI_asym 是新方向)

**实现**:
- GitHub: https://github.com/SBOZKAN/GNN
- 需要分子动力学预计算
- 实现复杂, ROI 较低

**风险**: 即使实现也可能因 OOD 失效, 重复 Round 2 教训

### 4.2 多任务学习 (2-3 天)

**目的**: 联合训练 brightness + Tm + fold_prob 三个预测器

**架构**:
- 共享 ESM-2 嵌入
- 3 个 XGBoost/MLP heads
- 数据: 14万条 CFPS + 文献 Tm + pLDDT

**价值**: Tm 预测能力 (替代 ThermoMPNN)

**风险**: 训练数据可能不够 (Tm 数据稀少)

### 4.3 进化搜索 (3-5 天)

**思路**: 用遗传算法/贝叶斯优化
- 起点: Round 4 v5 Top-6
- 变异: 随机替换 1-2 个突变为 ESM 推荐
- 适应度: pLDDT + chromophore cb + Tm 文献

**预期**: 100-200 代后可能找到 v5 没有的隐藏高分候选

### 4.4 多模态设计 (1 周+)

**整合**: 序列 + 结构 + MSA + 文献知识
- 工具: AlphaFold2 + ProteinMPNN + ESM3 + 序列进化
- 端到端学习: 输入 "我要高 Tm 高 Finit", 输出候选

**这是终极方向, 但 ROI 较低, 比赛时间内不建议**

---

## 五、风险评估与决策矩阵

### 5.1 风险 vs 收益决策表

| 决策 | 风险 | 收益 | 建议 |
|------|------|------|------|
| **改 v5 提交** | 破坏最优 | 边际提升 | ❌ 不要改 |
| **用 ThermoMPNN 重新排名** | 改变 Top-6 顺序 | 更准 | 🟡 如有时间可做 |
| **跑更多 MPNN** | 计算耗时 | 新高分 | 🟢 推荐 |
| **加 ESM3 候选** | API 风险 | 全新空间 | 🟡 谨慎 |
| **写 PDF** | 文档质量 | 比赛必交 | ✅ **必做** |
| **建 GitHub** | 时间 | 比赛必交 | ✅ **必做** |
| **联系 root** | 0 风险 | 高价值 | ✅ **必做** |

### 5.2 战略决策树

```
剩余时间?
├─ < 6h: 只做 PDF + GitHub
├─ 6-24h: PDF + GitHub + ThermoMPNN 验证
├─ 1-3d: PDF + GitHub + 200 MPNN + PROSS
├─ 3-7d: PDF + GitHub + 全套新工具
└─ > 7d: 重新做 Round 5 (新策略)
```

---

## 六、关键决策的"如果..."假设

### 如果比赛截止还有 1 天:
- ✅ 必做 PDF + GitHub
- ✅ 联系 root 询问细节
- 🟡 如果 v5 信心足, 不改提交
- ❌ 不做新工具部署

### 如果比赛截止还有 3 天:
- ✅ PDF + GitHub
- ✅ 跑 200 MPNN 候选
- ✅ 部署 ThermoMPNN
- ✅ 用 ThermoMPNN 重新排名
- 🟡 PROSS 在线服务 (1-2 小时)
- ❌ ESM3, ESM-IF1 (ROI 太低)

### 如果比赛截止还有 1 周:
- ✅ PDF + GitHub
- ✅ 跑 500 MPNN 候选 (5 个骨架, 5 温度)
- ✅ PROSS + ThermoMPNN + ESM-IF1 全套
- ✅ v6 重新排名
- 🟡 进化搜索 (3-5 天, 风险高)

### 如果发现 v5 Top-1 实际是 H1 (而非 MPNN):
- 重新跑 v5 选 H1 为 Seq 1
- MPNN_T01_014 降到 Seq 2 (作为最高 pLDDT 备份)
- 重新评估 v5 其它候选

### 如果 MPNN_T01_014 实际折叠失败 (pLDDT 误报):
- 这是低概率 (pLDDT 68 几乎确定能折叠)
- 但若发生, Top-6 中 5 条仍可正常工作
- 5 骨架多样性是保险

### 如果 CFPS 系统与文献差异大:
- mBaoJin (Tm 92°C) 实际可能 < 50°C (环境相关)
- sfGFP 风格 (Tm 80°C) 更可靠
- MPNN 候选 pLDDT 高但 Tm 估值不确定

---

## 七、给 Round 5 接手者的关键建议

### 7.1 三句真言

1. **"先交材料, 再追完美"** - PDF + GitHub 是硬门槛
2. **"pLDDT 是结构质量, 估值不等于真实"** - MPNN 高分有风险
3. **"比赛规则不明, 联系 root"** - 不确定就问

### 7.2 5 个必做

1. 联系 root 询问 CFPS 协议和截止日期
2. 写设计思路 PDF (含 LLM Agent 逻辑树)
3. 建立 GitHub 仓库 (含 6 份文档)
4. **保留 v5 提交** (不要轻易改)
5. 用新工具做 2-3 个候选 (高 ROI)

### 7.3 5 个不要

1. 不要用 ML 打分 OOD (Round 2 教训)
2. 不要在中文路径跑 ProteinMPNN
3. 不要追求 100% 完美 (BP Top-1 已足够)
4. 不要在 v5 的提交上做小修改 (破坏最优)
5. 不要在 Round 5 重新设计已有高分候选 (边际收益低)

### 7.4 Round 5 的"如果只有 1 个动作"是什么?

**如果只能做 1 件事**: **联系 root 询问比赛细节**
- 信息价值: 极高
- 时间成本: 5 分钟
- 决策影响: 可能完全改变 Round 5 策略

**如果能多做 1 件事**: **跑 200 MPNN 候选**
- 找到 2-3 条 pLDDT 65+ 的高分候选
- 替换 v5 的次优候选
- 显著降低比赛失败风险

**如果能多做 2 件事**: **部署 ThermoMPNN**
- 真实 ΔΔG 预测
- 替代间接 pLDDT 估值
- 更准的 Tm 排序

---

## 八、附录: 关键 URL 与联系方式

### 8.1 项目内

- 6 份文档: `docs/round4/round4_*.md`
- 最终提交: `work/round4/submission_round4_v5.csv`
- 候选池: `work/round4/final_6_round4_v5.json`
- MPNN 输出: `work/round4/mpnn_output_final/`

### 8.2 外部资源

| 资源 | URL |
|------|-----|
| ProteinMPNN | https://github.com/dauparas/ProteinMPNN |
| ThermoMPNN | https://github.com/Kuhlman-Lab/ThermoMPNN |
| ESM3 Forge | https://forge.evolutionaryscale.ai |
| ESM-IF1 | https://github.com/facebookresearch/esm |
| RCSB PDB | https://www.rcsb.org |
| PyTorch 国内镜像 | https://mirror.sjtu.edu.cn/pytorch-wheels/cu128 |
| htFuncLib 文献 | https://github.com/Fleishman-Lab |

### 8.3 关键时间点

| 时间 | 事件 |
|------|------|
| 2026-06-22 | Round 4 v5 完成 |
| ?-??-?? | 比赛截止 (**联系 root 确认**) |
| 截止前 24h | 最后冲刺 (PDF + GitHub + 不改提交) |
| 截止前 1h | 停止一切优化, 只做验证 |

---

**完成日期**: 2026-06-22
**作者**: Trae AI Agent
**最后更新**: Round 4 v5 提交就位
