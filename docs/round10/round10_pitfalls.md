# Round 10 踩坑记录

> **目标**: 总结 R10 期间遇到的真实问题及解决方案，避免后续 Agent 重蹈覆辙

---

## 1. ProteinMPNN 参数配置错误

### 问题
ProteinMPNN 文档说 `--fixed_positions_list` 指定固定残基，实际参数名是 `--fixed_positions_jsonl`，传入 `.list` 格式文件 MPNN 静默忽略。

### 后果
R4 实际为**无约束全自由设计**，与"5 核心固定"初衷不符。但反而效果最好（pLDDT=68.3）。

### 解决
- 严格按 ProteinMPNN `--help` 输出配置参数
- 验证 MPNN 输出文件的 `fixed_chains` 字段确认约束生效

### 教训
> **静默忽略比报错更危险**。每次 MPNN 跑完后必须核对输出文件。

---

## 2. 2B3P 晶体结构残缺

### 问题
sfGFP 晶体结构 2B3P 只有 228 个残基，**缺少 N 端 M1 和 C 端 6 个残基**。

### 后果
R10r2 用 2B3P + MPNN → 生成截短序列（231aa，N 端不是 M），全部不符合竞赛规则。

### 解决
改用 ESMFold 预测的 238aa 完整 PDB（虽然部分区域低置信，但完整）。

### 教训
> **骨架完整性 > 骨架精度**。MPNN 在不完整骨架上无法生成完整序列。

---

## 3. 固定 156 残基过度限制

### 问题
R10r1 固定所有 pLDDT ≥ 70 的 156 个残基 + 5 个核心，重设计 62 个低 pLDDT 残基。

### 后果
MPNN 设计空间过小，Top 1 排序分仅 0.679（vs R10r3 的 0.815）。

### 解决
改为仅固定 5 个核心残基（T65, Y66, G67, R96, E222），重设计 233 个残基。

### 教训
> **核心残基固定的最优数量在 5-10 之间**，固定 50+ 残基会大幅降低 MPNN 表达力。

---

## 4. ESMFold num_recycles 偏差

### 问题
R10 前期用 num_recycles=1 计算（快但不准），Top 1 score 0.83。改用 num_recycles=8 后，score 降至 0.815。

### 后果
低 recycles 下 ESMFold 给出虚高 pLDDT/pTM，**不严格符合竞赛规则**。

### 解决
所有最终评分用 num_recycles=8 重算。

### 教训
> **必须按竞赛规则 num_recycles=8 评估**，否则有"分数虚高"风险。

---

## 5. Windows 中文路径问题

### 问题
ProteinMPNN 在中文路径下 `torch.load` 崩溃。

### 解决
将 ProteinMPNN 复制到 `C:\proteinmpnn_r10`（纯英文路径）。

### 教训
> 关键工具链部署在英文路径，避免 Windows 中文路径陷阱。

---

## 6. ProteinMPNN 输出文件名含冒号

### 问题
MPNN 默认输出文件名含 `:`（如 `T=0.1`），Windows 不允许冒号。

### 解决
修改 `protein_mpnn_run.py` 中 `batch_clones[0]['name']` 使用 `os.path.basename()` 清理。

### 教训
> 跨平台路径需测试。Windows 工具部署到 Linux/Mac 之前应做完整路径回归测试。

---

## 7. C 端尾 pLDDT 难以提升

### 问题
R10 Top 6 的 200-238 残基 pLDDT 仍为 40-55，无论如何重设计都不提升。

### 根因
ESMFold 对 sfGFP 的 C 端尾（loop 区）预测本身置信度低 → 骨架不准 → MPNN 也学不准。

### 教训
> **预测 PDB 的低置信区域存在"骨架-序列"负反馈**。下一轮应对 C 端做截断或单独处理。

---

*踩坑记录作者: Trae AI Agent (Claude)*
*最后更新: 2026-06-23*
