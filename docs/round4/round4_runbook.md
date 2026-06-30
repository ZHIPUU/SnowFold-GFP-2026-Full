# Round 4 运行指南 (Runbook)

> **目的**: 提供从零复现 Round 4 全部工作的 step-by-step 步骤
> **适用**: 新接手 AI 完整复现 Round 4 / 准备 Round 5 改进
> **预计时间**: 完整复现 ~3-4 小时

---

## 零、复现前准备

### 0.1 检查环境 (5 分钟)

```powershell
# Python 版本
python --version  # 应为 3.14+

# PyTorch 是否 GPU 版
python -c "import torch; print(torch.__version__, torch.cuda.is_available())"
# 期望: 2.11.0+cu128 True
# 如果不是, 参考 pitfalls.md 坑 #1 修复

# GPU 是否识别
nvidia-smi  # 应看到 RTX 5080

# 项目目录存在
Test-Path "D:\生信\2026Protein Design\AAseqs of 5 GFP proteins_20260511.txt"
Test-Path "D:\生信\2026Protein Design\Exclusion_List.csv"
Test-Path "D:\生信\2026Protein Design\GFP_data.xlsx"
```

### 0.2 准备必要目录 (5 分钟)

```powershell
cd "D:\生信\2026Protein Design"
mkdir work\round4\pdbs
mkdir work\round4\logs
mkdir docs\round4
```

### 0.3 下载 5 个 GFP PDB (2 分钟)

```powershell
# 方式 1: 浏览器下载 https://files.rcsb.org/download/{ID}.pdb
# 方式 2: PowerShell
$pdbIds = @("2B3P","2WUR","7LG4","2HPW","8QBJ")
foreach ($id in $pdbIds) {
    $url = "https://files.rcsb.org/download/$id.pdb"
    $out = "D:\生信\2026Protein Design\work\round4\pdbs\$id.pdb"
    Invoke-WebRequest -Uri $url -OutFile $out -UseBasicParsing
    Write-Host "$id OK"
}
```

### 0.4 准备 ProteinMPNN (5 分钟)

```powershell
cd "D:\生信\2026Protein Design\work\round4"
git clone https://gitcode.com/gh_mirrors/pr/ProteinMPNN.git
# 如果 git 慢, 用国内镜像

# 创建纯英文临时目录并复制 MPNN
mkdir C:\Temp\mpnn_work
Copy-Item -Recurse ProteinMPNN C:\Temp\mpnn_work\
Copy-Item pdbs\2B3P.pdb C:\Temp\mpnn_work\
Copy-Item pdbs\2WUR.pdb C:\Temp\mpnn_work\
```

---

## 一、v1: 基础候选生成 (15 分钟)

### 1.1 准备 5 个 WT 序列 (1 分钟)

读取 `D:\生信\2026Protein Design\AAseqs of 5 GFP proteins_20260511.txt`, 提取:
- avGFP (238 aa)
- sfGFP (238 aa)
- amacGFP (238 aa)
- cgreGFP (235 aa)
- ppluGFP (222 aa)
- mBaoJin (234 aa, 来自 PDB 8QBJ 去掉 RS cloning artifact)

### 1.2 生成 v1 候选 (5 分钟)

```bash
cd "D:\生信\2026Protein Design"
python work/round4/01_seed_candidates.py
```

**预期输出**:
- `candidates_round4.json` (13 条)
- 全部通过验证 (220-250 aa, M 开头, chromophore 完整, 不在排除列表)

### 1.3 ESMFold 评估 v1 (5 分钟)

```bash
python work/round4/02_esmfold_screen.py
```

**预期输出**:
- `esmfold_round4.json` (13 条 pLDDT/pTM)
- 全部在 GPU 上跑 (~5s/条), 总耗时 ~70 秒

### 1.4 评分 v1 + 选 Top-6 (2 分钟)

```bash
python work/round4/03_score_and_select.py
```

**预期**: v1 Top-6 + 提交文件

---

## 二、v2: 多样性优化 (30 分钟)

### 2.1 多样性诊断 (2 分钟)

```bash
python work/round4/diagnose_diversity.py
```

**预期输出**: 6 候选的汉明距离矩阵
**关键发现**: Round 3 风格太相似, 5/6 汉明距离 1-4

### 2.2 生成多骨架候选 (10 分钟)

```bash
python work/round4/04_diversity_optimization.py
```

**预期**: 12 条新候选 (X1-X4 avGFP, Y1-Y3 cgreGFP, Z1-Z3 amacGFP, P1-P3 mBaoJin)

### 2.3 ESMFold 评估扩展池 (5 分钟)

```bash
python work/round4/05_esmfold_v2.py
```

**预期**: 33 条 pLDDT 数据

### 2.4 修正 mBaoJin 残基 (3 分钟)

发现 mBaoJin pos 173=K (不是 D), 173, 194, 222 等是 V/K/L (不是论文假设的 D/E/D)

```bash
python work/round4/01c_fix_mbaojin.py
```

**预期**: 8 条正确 mBaoJin 候选 (M1-M8)

### 2.5 v2 最终选择 + 预测 (3 分钟)

```bash
python work/round4/06_final_select_v2.py
python work/round4/07_score_estimate.py
```

**预期**: 6 候选覆盖 4 骨架 (sfGFP×2, avGFP×2, amacGFP×1, mBaoJin×1)

---

## 三、v3: ProteinMPNN 引入 (60 分钟)

### 3.1 部署 ProteinMPNN (10 分钟)

```bash
# 在 C:\Temp\mpnn_work 中运行 ProteinMPNN
cd C:\Temp\mpnn_work
python ProteinMPNN/helper_scripts/parse_multiple_chains.py \
    --input_path . --output_path parsed.jsonl

# 检查解析结果
python -c "
import json
with open('parsed.jsonl') as f:
    p = json.loads(f.readline())
print('seq:', p['seq_chain_A'][:50])
print('len:', len(p['seq_chain_A']))
"
```

**预期**: 解析出 sfGFP 2B3P chain A 序列 (231 aa, 缺 N 端 M)

### 3.2 跑 sfGFP MPNN (3 个温度) (20 分钟)

```bash
cd C:\Temp\mpnn_work
# 写并运行完整脚本 (见 09e_proteinmpnn_v5.py)
python D:\生信\2026Protein Design\work\round4\09e_proteinmpnn_v5.py
```

**预期输出**:
- T01: 25 序列
- T03: 25 序列
- T05: 17 序列
- 总 67 条候选

**注意**: 第一次跑可能遇到"中文路径"或"输出路径不识别"错误, 参考 pitfalls.md 坑 #2 修复

### 3.3 加载 MPNN 候选 + 填回 X 占位符 (5 分钟)

```bash
python D:\生信\2026Protein Design\work/round4/10c_load_mpnn_v3.py
```

**预期**: 12 条 MPNN 候选保存到 `mpnn_candidates_final.json`

**关键**: 正确处理 X 占位符, 补回 N 端 M, 验证 chromophore 存在

### 3.4 ESMFold 评估 MPNN 候选 (3 分钟)

```bash
python D:\生信\2026Protein Design\work/round4/11_esmfold_mpnn.py
```

**预期输出**: 12 条 pLDDT 数据
**惊喜**: MPNN_T01_014 pLDDT 68.3 = 全场最高!

### 3.5 v3 最终选择 (3 分钟)

```bash
python D:\生信\2026Protein Design\work/round4/12_final_select_v3.py
```

**预期**: 6 候选 (5 骨架: sfGFP, avGFP, sfGFP_MPNN, amacGFP, mBaoJin)

---

## 四、v4: avGFP MPNN (45 分钟)

### 4.1 跑 avGFP MPNN (15 分钟)

```bash
cd C:\Temp\mpnn_work
# 复制 avGFP PDB 到工作目录
Copy-Item D:\生信\2026Protein Design\work\round4\pdbs\2WUR.pdb C:\Temp\mpnn_work\avGFP.pdb

# 写并运行 (见 14_mpnn_multi_scaffold.py, 15_fix_avgfp_mpnn.py)
python D:\生信\2026Protein Design\work/round4/14_mpnn_multi_scaffold.py
python D:\生信\2026Protein Design\work/round4/15_fix_avgfp_mpnn.py
```

**注意**:
- 2WUR chromophore TYG 在 1-based pos 36 (不是 65)
- 必须用解析后真实位置, 不能套用 sfGFP 编号
- 详见 pitfalls.md 坑 #5

### 4.2 加载 avGFP MPNN 候选 (3 分钟)

```bash
python D:\生信\2026Protein Design\work/round4/17_load_avgfp_mpnn.py
```

**预期**: 8 条 avGFP MPNN 候选

### 4.3 ESMFold 评估 avGFP MPNN (2 分钟)

```bash
python D:\生信\2026Protein Design\work/round4/18_esmfold_av_mpnn.py
```

**预期输出**: 8 条 pLDDT
**亮点**: MPNN_av_T03_v2_001 pLDDT 61.5

### 4.4 v4 最终选择 (3 分钟)

```bash
python D:\生信\2026Protein Design\work/round4/19_final_select_v4.py
```

**预期**: 6 候选覆盖 6 骨架 (新增 avGFP_MPNN)

---

## 五、v5: 修正 Tm 估值 (10 分钟)

### 5.1 v5 最终选择 (5 分钟)

```bash
python D:\生信\2026Protein Design\work/round4/21_final_v5_with_corrected_tm.py
```

**预期输出**:
- `submission_round4_v5.csv` ⭐ **最终推荐提交**
- `final_6_round4_v5.json`

**关键变化**:
- MPNN 候选 Tm 估值按 pLDDT 分级 (88/84/81/78/72)
- 综合分排序重排

### 5.2 验证最终提交 (1 分钟)

```powershell
Get-Content "D:\生信\2026Protein Design\work\round4\submission_round4_v5.csv"
```

**验证**:
- 6 行数据
- 每行 Seq_ID 1-6
- 每条序列 234-238 aa
- 全部 M 开头
- 全部含 TYG 或 GYG

---

## 六、生成设计思路 PDF (1-2 小时)

### 6.1 收集内容 (15 分钟)

从以下文档/文件提取:
- `docs/round4/round4_report.md` (主报告)
- `docs/round4/round4_design_rationale.md` (设计思路)
- `docs/round4/round4_tech_reference.md` (技术参考)
- `work/round4/final_6_round4_v5.json` (最终候选)
- `work/round4/submission_round4_v5.csv` (提交 CSV)

### 6.2 写 PDF 生成脚本 (30 分钟)

参考 Round 1 的 `work/phase5/md_to_pdf.py`, 用 reportlab 生成。

```python
# 伪代码
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors

doc = SimpleDocTemplate("design_doc_round4.pdf", pagesize=A4)
styles = getSampleStyleSheet()
story = []

# 封面
story.append(Paragraph("Round 4 GFP Design Document", styles["Title"]))
story.append(Paragraph("YourTeamName", styles["Normal"]))
story.append(Paragraph("2026-06-22", styles["Normal"]))
story.append(Spacer(1, 0.5*inch))

# 1. 任务理解
story.append(Paragraph("1. Task Understanding", styles["Heading1"]))
story.append(Paragraph("目标: 设计兼具高荧光亮度和优良热稳定性的 GFP 变体", styles["Normal"]))
story.append(Paragraph("综合分 = Finit/Finit_WT × Ffinal/Finit (72°C 热处理)", styles["Normal"]))
# ...

# 2. 方法论
story.append(Paragraph("2. Methodology", styles["Heading1"]))
# ...

# 3. 6 候选详解
story.append(Paragraph("3. Final 6 Candidates", styles["Heading1"]))
# 表格 + 序列

doc.build(story)
```

### 6.3 必含内容 (1 小时写)

按比赛要求 5 大块:

1. **任务理解与目标拆解** (10%)
2. **算法管线** (25%): 数据→加性模型→ESMFold→ProteinMPNN→评分
3. **如何平衡亮度和稳定性目标** (15%)
4. **筛选 6 候选的原因** (20%)
5. **LLM Agent 逻辑树 + 关键日志** (30%): 完整复现路径

---

## 七、建立 GitHub 开源仓库 (30 分钟)

### 7.1 创建仓库 (5 分钟)

1. 在 https://github.com 创建仓库 `2026-protein-design`
2. README.md 写项目说明
3. 关联本地仓库

### 7.2 准备仓库内容 (15 分钟)

```bash
cd D:\生信\2026Protein Design
# 创建 .gitignore
echo "__pycache__/" > .gitignore
echo "work/round4/mpnn_output*/" >> .gitignore
echo "work/round4/ProteinMPNN/" >> .gitignore  # 太大, 不上传
echo "C:\Temp/" >> .gitignore

# 重要文件
git add README.md
git add docs/round4/  # 6 份文档
git add work/round4/01_seed_candidates.py  # 关键脚本
git add work/round4/06_final_select_v2.py
git add work/round4/12_final_select_v3.py
git add work/round4/19_final_select_v4.py
git add work/round4/21_final_v5_with_corrected_tm.py
git add work/round4/submission_round4_v5.csv  # 最终提交
git add requirements.txt

git commit -m "Round 4 v5 final submission"
git push
```

### 7.3 README.md 必含 (10 分钟)

```markdown
# 2026 Protein Design — GFP Variant Engineering

## 项目背景
SynBio Challenges 2026 — Protein Design by AI Track

## 环境配置
- Python 3.14
- PyTorch 2.11.0+cu128
- RTX 5080 (16GB VRAM)
- 详见 docs/round4/round4_tech_reference.md

## 快速运行
```bash
# 安装依赖
pip install torch torchvision --index-url https://mirror.sjtu.edu.cn/pytorch-wheels/cu128
pip install transformers pandas numpy

# 下载 PDB (5 个)
# 见 docs/round4/round4_runbook.md 0.3 节

# 复现 Round 4
python work/round4/01_seed_candidates.py
python work/round4/02_esmfold_screen.py
# ... 详见 runbook

# 提交文件: work/round4/submission_round4_v5.csv
```

## 详细文档
- docs/round4/round4_report.md (主报告)
- docs/round4/round4_design_rationale.md (设计思路)
- docs/round4/round4_tech_reference.md (技术参考)
- docs/round4/round4_pitfalls.md (踩坑指南)
- docs/round4/round4_runbook.md (运行指南)
- docs/round4/round4_next_steps.md (下一步)

## 关键结果
- 6 候选 (sfGFP/avGFP/sfGFP_MPNN/avGFP_MPNN/amacGFP/mBaoJin)
- 最高 pLDDT 68.3 (MPNN_T01_014)
- 最高 pTM 0.765 (MPNN_T01_014)
- 预测综合分 1.23 (H1_avGFP_sfGFP_acid3_I152S)
```

---

## 八、完整流程时间表

| 阶段 | 步骤 | 时间 | 总累计 |
|------|------|------|--------|
| 准备 | 环境检查 | 5min | 5min |
| 准备 | 下载 PDB | 2min | 7min |
| 准备 | 部署 ProteinMPNN | 5min | 12min |
| v1 | 候选生成 | 5min | 17min |
| v1 | ESMFold 评估 | 5min | 22min |
| v1 | 评分选 Top-6 | 2min | 24min |
| v2 | 多样性诊断 | 2min | 26min |
| v2 | 多骨架候选 | 10min | 36min |
| v2 | ESMFold 评估 | 5min | 41min |
| v2 | mBaoJin 修正 | 3min | 44min |
| v2 | v2 选择 + 预测 | 3min | 47min |
| v3 | 跑 sfGFP MPNN | 20min | 67min |
| v3 | 加载 MPNN | 5min | 72min |
| v3 | ESMFold MPNN | 3min | 75min |
| v3 | v3 选择 | 3min | 78min |
| v4 | 跑 avGFP MPNN | 15min | 93min |
| v4 | 加载 avGFP MPNN | 3min | 96min |
| v4 | ESMFold avGFP MPNN | 2min | 98min |
| v4 | v4 选择 | 3min | 101min |
| v5 | v5 选择 + 预测 | 5min | 106min |
| 提交 | PDF 生成 | 60min | 166min |
| 提交 | GitHub 仓库 | 30min | 196min |
| **总** | **完整复现 + 提交** | **~3.3 小时** | |

---

## 九、验证清单 (复现完成前对照)

- [ ] `work/round4/submission_round4_v5.csv` 存在且符合格式
- [ ] 6 候选长度 234-238 aa
- [ ] 全部 M 开头
- [ ] 全部含 TYG/SYG/GYG
- [ ] 全部不在 Exclusion_List 中
- [ ] 6 候选来自 ≥ 5 个不同骨架 (多样性)
- [ ] docs/round4/ 6 份文档齐全
- [ ] GitHub 仓库已推送
- [ ] README.md 含完整环境配置和运行说明
- [ ] 设计思路 PDF 已生成
- [ ] (可选) ThermoMPNN Tm 验证
- [ ] (可选) PROSS 在线服务

---

## 十、紧急恢复 (如果出错)

### 10.1 找不到最佳提交

```bash
# 列出所有 v* 提交
ls "D:\生信\2026Protein Design\work\round4\submission_round4_v*.csv"

# 推荐顺序:
# v5 (修正 Tm) > v4 (含 avGFP_MPNN) > v3 (含 sfGFP_MPNN) > v2 (多样性)
```

### 10.2 重新跑 ESMFold 但 GPU 失败

```python
# 强制 CPU
import os
os.environ["CUDA_VISIBLE_DEVICES"] = ""  # 禁用 GPU
# 加载模型时会自动用 CPU
# 每序列 130s (慢 26x 但能跑)
```

### 10.3 MPNN 文件丢失

```bash
# 重新克隆
cd C:\Temp
rm -rf mpnn_work
mkdir mpnn_work
cd mpnn_work
git clone https://gitcode.com/gh_mirrors/pr/ProteinMPNN.git
# 复制 PDB 等
```

### 10.4 排除列表太大加载慢

```python
# 用更高效的方法
import pandas as pd
excl = pd.read_csv("Exclusion_List.csv", usecols=["Sequence"], dtype=str)
excl_seqs = set(excl["Sequence"].str.strip())
# 单列读入节省内存
```

---

**完成日期**: 2026-06-22
**作者**: Trae AI Agent
**总时间**: ~3-4 小时 (从零复现)
