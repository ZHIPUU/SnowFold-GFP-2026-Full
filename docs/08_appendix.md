# 08 · 附录(Appendix)

> 命令速查、文件清单、环境配置、模型训练命令、引用。

---

## 一、文件清单

### 根目录
| 文件 | 用途 |
|---|---|
| `README.md` | 项目入口 |
| `submission_yourteamname.csv` | **当前 Round 2 提交**(6 条候选) |
| `AAseqs of 5 GFP proteins_20260511.txt` | 5 个 GFP WT 序列 |
| `Exclusion_List.csv` | 13.5 万条排除列表 |
| `GFP_data.xlsx` | 官方 14 万条 CFPS 亮度数据 |
| `submission_template.csv` | 提交模板 |
| `requirements.txt` | Python 依赖 |
| `run_all.py` | Round 1 一键复现 |
| `2026Protein Design in Synbio challenges.pdf` | 比赛规则 |
| `Basic Tutorial on Protein Design.ipynb` | 比赛给的 tutorial |

### docs/(本文档)
| 文件 | 主题 |
|---|---|
| `01_achievements.md` | 两轮成果汇总 |
| `02_methodology.md` | 完整 pipeline |
| `03_challenges.md` | 难点与坑 |
| `04_open_questions.md` | 待解疑点 |
| `05_next_steps.md` | 下一步方向 |
| `06_paper_kb.md` | 论文突变知识库 |
| `07_handoff.md` | 接手指南 |
| `08_appendix.md` | 本附录 |

### work/round2/(Round 2 主目录)
| 文件 | 用途 |
|---|---|
| `step1_esm650m_embeddings.py` | Step 1: ESM2-650M 嵌入 141K |
| `step1_log.txt` | Step 1 日志 |
| `esm650m_embeddings.npy` | 141K × 1280 嵌入(0.7 GB) |
| `esm650m_ids.csv` | 嵌入对应的元数据 |
| `all_variants.pkl` | 141K 变体的完整序列(pickle) |
| `all_variants.csv` | 同上,CSV 格式 |
| `build_variants.py` | 从突变字符串重建完整序列 |
| `step2_train.py` | Step 2: XGBoost GPU 训练 |
| `step2_xgboost_gpu.model` | Step 2 模型(27 MB) |
| `step2_xgboost_gpu.json` | Step 2 模型指标 |
| `step2_baseline_ridge.json` | Ridge baseline 指标 |
| `step2_log.txt` | Step 2 日志 |
| `step3a_score_141k.py` | Step 3a: 用 XGBoost 打分 141K |
| `step3a_141k_scored.csv` | 141K 评分结果 |
| `step3a_summary.json` | Step 3a 摘要 |
| `step3a_finit_rel.py` | Finit/Finit_WT 计算 |
| `step3b_embed_top10k.py` | Step 3b: Top 10K 重嵌入(已废弃) |
| `step3b_summary.json` | Step 3b 摘要(显示全 0.000×) |
| `step3c_design_candidates.py` | Step 3c: 11 条手工候选设计 |
| `candidates_round2_design.csv` | **11 条候选** |
| `step4_esmfold.py` | Step 4: ESMFold(已废弃,SSL 卡住) |
| `step4_validate_no_esmfold.py` | Step 4 替代:轻量级验证 |
| `step4_validation.csv` | 候选评分 |
| `step5_generate_submission.py` | Step 5: 生成 paper-knowledge Top-6 提交 |
| `step5b_revert_submission.py` | Step 5b: 回退到 paper-knowledge 提交 |
| `submission_round2.csv` | Step 5 生成的中间提交 |
| `stepC_epistasis.py` | **Step C: Epistasis 模型(R²=0.92)** |
| `stepC_xgboost_epistasis.model` | **Epistasis 模型**(2.8 MB) |
| `stepC_summary.json` | Step C 模型指标 |
| `stepC2_score_candidates.py` | Step C2: 候选打分(显示 OOD 失效) |
| `stepC2_candidates_scored.csv` | Step C2 评分结果 |
| `top10k_*` | Step 3b 中间产物(已废弃) |
| `step*_stdout.log`, `step*_stderr.log` | 各步骤运行日志 |

### work/phase1-6/(Round 1 各阶段)
- phase1:数据探索 + 加性 Ridge
- phase2:ESM2-150M PLL 评分
- phase3:终选 6 条
- phase5:设计思路 PDF
- phase6:综合分预测

### referencepaper/(5 篇 PDF)
- superfolder.pdf (628 KB)
- TGP.pdf (2.4 MB)
- StayGold.pdf (3.2 MB)
- mBaoJin.pdf (1.4 MB)
- nature-Local fitness landscape.pdf (3.1 MB)

---

## 二、环境配置

### 硬件
- GPU:**RTX 5080 Laptop**(sm_120, Blackwell, 16 GB VRAM)
- 单卡,不并发

### Python 环境
```bash
# Python 3.14
# PyTorch nightly + CUDA 12.8
pip install --pre torch --index-url https://download.pytorch.org/whl/nightly/cu128

# 关键库
pip install fair-esm xgboost optuna lightgbm scikit-learn pandas numpy scipy openpyxl transformers
```

### 已知问题
- LightGBM 4.6.0 **无 cu128 预编译**,只能 CPU
- XGBoost 3.3.0 **有 cu128 预编译**,GPU OK
- fair-esm 2.0 **不含 ESMFold**(需要 esmfold_v1 单独下载)
- HuggingFace `huggingface_hub` 在 Python 3.14 上 SSL 卡住

---

## 三、关键命令速查

### 重新嵌入(Step 1, ~83 min)
```bash
python "D:\生信\2026Protein Design\work\round2\step1_esm650m_embeddings.py"
```

### 重新训练 XGBoost(Step 2, ~2 min)
```bash
python "D:\生信\2026Protein Design\work\round2\step2_train.py"
```

### 重新训练 Epistasis 模型(Step C, ~3 min)
```bash
python "D:\生信\2026Protein Design\work\round2\stepC_epistasis.py"
```

### 重新打分 11 候选(Step C2, ~10 sec)
```bash
python "D:\生信\2026Protein Design\work\round2\stepC2_score_candidates.py"
```

### 重新生成 11 候选(Step 3c, ~1 sec)
```bash
python "D:\生信\2026Protein Design\work\round2\step3c_design_candidates.py"
```

### 重新生成提交(Step 5b, ~1 sec)
```bash
python "D:\生信\2026Protein Design\work\round2\step5b_revert_submission.py"
```

### 验证嵌入可加载
```python
import numpy as np
emb = np.load(r"D:\生信\2026Protein Design\work\round2\esm650m_embeddings.npy", mmap_mode="r")
print(emb.shape, emb.dtype)  # (141144, 1280) float32
```

### 验证提交合规
```python
import pandas as pd
sub = pd.read_csv(r"D:\生信\2026Protein Design\submission_yourteamname.csv")
for _, r in sub.iterrows():
    s = r["Sequence"]
    print(f"Seq {r['Seq_ID']}: len={len(s)} M-start={s[0]=='M'} chromo={'✓' if any(t in s for t in ['TYG','SYG','GYG']) else '✗'}")
```

### 全量排除列表检查
```python
import pandas as pd
excl = pd.read_csv(r"D:\生信\2026Protein Design\Exclusion_List.csv")
excl_seqs = set(excl.iloc[:, 0].astype(str))
sub = pd.read_csv(r"D:\生信\2026Protein Design\submission_yourteamname.csv")
for _, r in sub.iterrows():
    print(f"Seq {r['Seq_ID']}: {'❌' if r['Sequence'] in excl_seqs else '✓'}")
```

---

## 四、模型训练命令详解

### Step 2 (XGBoost GPU)
```python
params = {
    "objective": "reg:squarederror",
    "tree_method": "hist",
    "device": "cuda",
    "max_depth": 8,
    "learning_rate": 0.05,
    "subsample": 0.85,
    "colsample_bytree": 0.7,
    "min_child_weight": 5,
    "reg_alpha": 0.1,
    "reg_lambda": 1.0,
}
# early_stop=50, num_boost_round=2000
# best_iter=1775, 训练 84s
```

### Step C (Epistasis + XGBoost GPU)
- 同上,但增加 9 维特征:
  - type one-hot (4)
  - n_mut (1)
  - single_sum, single_mean, single_max, single_min (4)
- best_iter=193, 训练 ~3 min

---

## 五、WT 序列速查(从 `AAseqs of 5 GFP proteins_20260511.txt`)

| 名称 | 长度 | 关键特征 |
|---|---|---|
| sfGFP | 238 | 含 F64L/S65T/F99S/M153T/V163A 等 11 突变,Tm ~78°C |
| avGFP | 238 | GFP 经典,Tm ~64°C,亮度 baseline 3.72 log10 |
| amacGFP | 238 | A. macrodactyla GFP,Tm ~?,亮度 baseline 3.97 log10 |
| cgreGFP | 235 | Cytaeis gregaria GFP,**亮度 baseline 4.50 log10 (最高)** |
| ppluGFP | 222 | Pontellina plumata GFP,Tm ~78°C,亮度 baseline 4.23 log10 |

### WT 关键位置(用户早期纠正后)
- avGFP: pos 65=S, 152=I, 163=V, 171=I, 206=A
- cgreGFP: pos 65=S, 130=G, 167=N, 192=W
- ppluGFP: pos 137=T, 159=S, 171=D, 199=L
- amacGFP: pos 166=K(经典 I166T hotspot)

### WT baseline brightness (log10)
- avGFP: 3.72
- amacGFP: 3.97
- cgreGFP: 4.50
- ppluGFP: 4.23
- sfGFP: ~4.17(估计)

---

## 六、引用文献

### 数据集
- Sarkisyan KS, Bolotin DA, Meer MV, et al. Local fitness landscape of the green fluorescent protein. Nature. 2016;533(7603):397-401.

### 参考论文
1. Pédelacq JD, Cabantous S, Tran T, Terwilliger TC, Waldo GS. Engineering and characterization of a superfolder green fluorescent protein. Nat Biotechnol. 2006;24(1):79-88.
2. Close DW, Don Paul C, Langan PS, et al. TGP, an extremely stable, non-aggregating fluorescent protein created by structure-guided surface engineering. Proteins. 2015;83(7):1225-1237.
3. Hirano M, Ando R, Shimozono S, et al. A highly photostable and bright green fluorescent protein. Nat Biotechnol. 2022;40(8):1132-1142.
4. (mBaoJin 论文,待精确引用)
5. (Local fitness landscape 已引)

### 模型架构
- ESM2: Lin Z, et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. Science. 2023;379(6637):1123-1130.
- XGBoost: Chen T, Guestrin C. XGBoost: A Scalable Tree Boosting System. KDD '16.

### 嵌入文件位置
- ESM2-650M 模型:`C:\Users\A\.cache\torch\hub\checkpoints\esm2_t33_650M_UR50D.pt.75c71e41769c4391ba0186bc8c92d0f7.partial`
- ESM2-650M 嵌入(本地):`D:\生信\2026Protein Design\work\round2\esm650m_embeddings.npy`

---

## 七、版本控制

- 项目未使用 git(直接放在用户工作目录)
- 大文件(嵌入、模型)未 commit,文档化位置
- `submission_yourteamname.csv` 是当前最终交付

---

## 八、已知 Bug 与注意事项

1. **路径中文乱码**:PowerShell 调 `Get-Content` 中文路径经常乱码,用 `-LiteralPath` + `-Encoding UTF8`
2. **Start-Process 路径截断**:用 `subprocess.Popen` 替代
3. **HF SSL 卡**:用 `urllib` 替代 `huggingface_hub`
4. **Embedding 加载 ~570 missing keys**:可忽略(strict=False)
5. **PowerShell `$_` 在 bash 工具被解析**:用 Python psutil 替代
6. **RTX 5080 sm_120**:必须用 PyTorch nightly cu128

---

## 九、联系信息

- **AI session** (mvs_28b21be1eb074dfa9bcbf4c28733a5fa) - 当前 session
- **Root session** (mvs_7706ec6fc7c74872a207498ef6d551ed) - 用户
- **比赛**:Synbio Challenges 2026 — Protein Design

---

*如果有任何不清楚的地方,先看 `docs/03_challenges.md` 的对应章节。*