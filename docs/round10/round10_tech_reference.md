# Round 10 技术参考

> **目标**: 提供 R10 所需的关键技术参数、API、性能数据

---

## 一、核心工具版本

| 工具 | 版本 | 来源 | 用途 |
|:-----|:----:|:-----|:----|
| ProteinMPNN | v_48_020 | Baker Lab (GitHub) | 逆折叠 |
| ESMFold | facebook/esmfold_v1 | Meta AI (HuggingFace) | 结构预测 |
| PyTorch | 2.11.0+cu128 | pip | 深度学习框架 |
| transformers | 4.48.1 | HuggingFace | 模型加载 |
| Python | 3.14.0 | python.org | 运行时 |

## 二、硬件环境

```
GPU: NVIDIA GeForce RTX 5080 Laptop GPU
VRAM: 16GB
CUDA: 12.8
Driver: 570.x
RAM: 32GB
OS: Windows 11
```

## 三、性能指标

| 模型 | 输入 | 单次耗时 | VRAM | 批大小 |
|:-----|:-----|:--------:|:----:|:------:|
| ProteinMPNN | 238aa PDB | ~5s/50条 | 2GB | 25 |
| ESMFold (r=1) | 239aa | ~5s | 3GB | 1 |
| ESMFold (r=8) | 239aa | ~15s | 4GB | 1 |
| ESMFold (r=20) | 239aa | ~40s | 5GB | 1 |

## 四、关键参数

### 4.1 ProteinMPNN 参数

```bash
--pdb_path <path>           # 输入 PDB
--out_folder <path>         # 输出目录
--num_seq_per_target 50     # 每温度候选数
--sampling_temp "0.1 0.2 0.3 0.5 0.7"  # 多温度
--seed 37                   # 随机种子
--batch_size 25             # 批大小
--fixed_positions_jsonl <file>  # 固定残基（注意是 jsonl）
```

### 4.2 ESMFold 参数

```python
model = EsmForProteinFolding.from_pretrained(
    "facebook/esmfold_v1",
    low_cpu_mem_usage=True,
    local_files_only=True
).cuda()
model.trunk.set_chunk_size(128)  # 受显存限制
model.eval()

# num_recycles 关键
with torch.no_grad():
    output = model(tokens, num_recycles=8)  # 竞赛规则要求
```

### 4.3 排序分计算

```python
score = (
    0.40 * ptm +                          # pTM 权重 40%
    0.30 * (plddt_mean / 100) +           # 全局 pLDDT 权重 30%
    0.30 * (plddt_chromo / 100)           # 生色团 pLDDT 权重 30%
)
```

## 五、关键脚本

### 5.1 r10_iterative_repair.py — 主修复管线

**职责**: 接收 R4 MPNN_T01_014 序列，运行完整迭代修复流程

**关键步骤**:
1. ESMFold 预测 R4 序列结构
2. 残基级 pLDDT 分析
3. 生成 5 核心残基固定文件
4. 调用 ProteinMPNN 多温度采样
5. ESMFold num_recycles=8 验证所有候选
6. 输出 Top 20 JSON

### 5.2 r10_step1_analyze_best.py — R4 残基分析

**职责**: 对 R4 MPNN_T01_014 做残基级 pLDDT 分布分析

**输出**: 每个残基的 pLDDT 值，标识低置信区域

## 六、数据格式

### 6.1 候选 JSON 字段

```json
{
  "name": "R10r3_T02_078",
  "temperature": 0.2,
  "seq": "MPIPGDELLS...",
  "length": 239,
  "n_muts": 109,
  "plddt_mean": 76.69,
  "plddt_chromo_region": 81.57,
  "ptm": 0.8467,
  "sort_score": 0.8157,
  "fold_time_s": 5.5,
  "pass_ptm": true,
  "pass_plddt": true,
  "pass_chromo": true,
  "all_pass": true
}
```

### 6.2 提交 CSV 格式

```csv
Team_Name,Seq_ID,Sequence
SnowFold,1,MPIPGDELLS...
SnowFold,2,MAIPGDELLK...
...
```

## 七、关键文件路径

```
D:\生信\2026Protein Design\
├── work\round10\
│   ├── mpnn_output\                  # MPNN 输出
│   │   ├── parsed_chains.jsonl
│   │   └── fixed_positions.jsonl
│   ├── r10_iterative_repair.py
│   ├── r10_step1_analyze_best.py
│   ├── r10_round3_all.json           # 285 候选
│   ├── r10_round3_top20.json
│   ├── r10_round2_all.json           # R10r2 失败记录
│   ├── r10_round2_top30.json
│   ├── r10_top20.json
│   ├── r10_all_candidates.json
│   ├── t01_014_analysis.json         # R4 残基分析
│   ├── final_6_r10.json              # 最终 6 条
│   └── submission_r10_final.csv
├── C:\proteinmpnn_r10\              # MPNN 部署
└── 预选序列\
    ├── final_top6.json               # R10 最终
    ├── submission_top6.csv
    └── 实验总结报告.md
```

## 八、常见错误码

| 错误 | 原因 | 解决 |
|:-----|:----|:-----|
| `UnicodeDecodeError` | 中文路径 | 改用英文路径 |
| `RuntimeError: CUDA OOM` | batch_size 太大 | 减小至 1 |
| `KeyError: 'plddt'` | ESMFold 输出结构变化 | 检查 transformers 版本 |
| `np.int` 已移除 | numpy 2.x | 升级或替换代码 |

## 九、参考资料

- ProteinMPNN: https://github.com/dauparas/ProteinMPNN
- ESMFold: https://huggingface.co/facebook/esmfold_v1
- 2B3P PDB: https://www.rcsb.org/structure/2B3P
- 竞赛规则: `D:\生信\2026Protein Design\.trae\rules\目标及规范.md`

---

*技术参考作者: Trae AI Agent (Claude)*
*最后更新: 2026-06-23*
