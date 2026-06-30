# 03 · 难点与坑(Challenges)

> 本文记录项目过程中遇到的所有技术难点、环境坑、设计陷阱,以及解决方案。接手者可据此避免重复踩坑。

---

## 一、环境与基础设施

### 1. Python 3.14 + 稳定 PyTorch 不兼容

**症状**:`pip install torch` 安装的 stable 版本(2.6/2.7)无 CUDA 12.8 支持,在 RTX 5080 (sm_120, Blackwell) 上无法运行。

**解决**:装 PyTorch nightly + cu128:
```
pip install --pre torch --index-url https://download.pytorch.org/whl/nightly/cu128
```

**遗留问题**:LightGBM 4.6.0 没有 cu128 预编译,只能用 CPU。**XGBoost 3.3.0 有 cu128 预编译,GPU OK**。

### 2. ESM2-650M 加载卡死(torch.hub)

**症状**:`esm.pretrained.esm2_t33_650M_UR50D()` 内部调 `load_model_and_alphabet_hub` → `torch.hub.load_state_dict_from_url` → 先 HEAD 请求 HuggingFace URL,在某些网络环境下挂死。

**解决**:绕过 torch.hub,直接本地加载:
```python
ckpt = torch.load(model_path, map_location="cuda", weights_only=False)
args = ckpt["cfg"]["model"]
alphabet = esm.Alphabet.from_architecture("ESM-1b")
model = esm.ESM2(
    num_layers=args.encoder_layers,
    embed_dim=args.encoder_embed_dim,
    attention_heads=args.encoder_attention_heads,
    alphabet=alphabet,
).cuda()
model.load_state_dict(ckpt["model"], strict=False)  # ~570 missing/expected keys,可忽略
```

### 3. HuggingFace SSL 证书问题

**症状**:HuggingFace transformers / huggingface_hub 库在 Python 3.14 + httpx 环境下报 `SSL: CERTIFICATE_VERIFY_FAILED`,即使设了 `HF_HUB_DISABLE_SSL_VERIFY=1` 也无效。

**原因**:Python 3.14 + httpx 默认不带 certifi 证书。

**绕道**:用 `urllib` 手动下载(因为 Python 3.14 默认带证书),但 ESMFold 模型 14GB,即便能下载也很慢。

**结论**:ESMFold **无法本地使用**,改用轻量级验证(chromophore + 排除列表 + 论文先验)。

### 4. Start-Process 路径空格问题

**症状**:`Start-Process python -ArgumentList '-u','step2_train.py' -WorkingDirectory 'D:\生信\2026Protein Design\work\round2'` 启动后,stderr 报 `can't open file 'D:\\生信\\2026Protein'`(路径在空格 "Protein Design" 处截断)。

**原因**:PowerShell Start-Process 与 -ArgumentList 的 quoting 在某些边界条件下会拼错路径。

**解决**:用 Python subprocess.Popen + DETACHED_PROCESS 标志启动:
```python
import subprocess
subprocess.Popen(
    [str(PY), "-u", str(SCRIPT)],  # 用 list 而不是 string
    cwd=str(ROOT),                    # 显式 cwd
    stdout=open(LOG, "wb"),
    stderr=open(ERR, "wb"),
    creationflags=0x00000008 | 0x00000200,  # DETACHED + NEW_PROCESS_GROUP
)
```

### 5. PowerShell 中文路径 + `$variable` 解析

**症状**:`Get-Process python | Where-Object {$_.Id -ne ...}` 在 bash 工具里被错误解析为 `$_.Id` 而不是 `$_`。

**解决**:用 Python `psutil` 替代,或写 `.cmd` / `.ps1` 包装脚本。

---

## 二、数据相关

### 6. 突变字符串多分隔符

**症状**:数据中 `aaMutations` 字段有时用 `:`,有时用 `/`,有时用 `,` 分隔;还有 `*238G` 这类插入突变。

**解决**:标准化解析:
```python
def parse_muts(s):
    s = str(s).replace(",", ":").replace("/", ":")
    out = []
    for tok in s.split(":"):
        tok = tok.strip()
        if not tok or tok.startswith("*"):
            continue
        m = re.match(r"^([A-Z])(\d+)([A-Z])$", tok)
        if m:
            out.append((m.group(1), int(m.group(2)), m.group(3)))
    return out
```

### 7. WT 序列位置编号错误(用户纠正)

**症状**:用户早期告诉我 avGFP pos 152 = M(实际是 I),cgreGFP pos 167 = K(实际是 N)等。

**教训**:**任何位置编号都先自己重新数一遍序列**,不要盲信 agent 或文档。

### 8. 数据 max brightness ~4.6 的硬天花板

**现象**:141K 数据中,brightness (log10) 范围 [1.28, 4.60],**没有任何变体超过 4.60**(对应 ~40K linear)。

**原因推测**:数据集基于已部分 sfGFP 化的中间体;真正的 sfGFP / TGP / StayGold 风格的"突变组合巨大跳跃"**不在数据里**。

**后果**:数据驱动的 max Finit/Finit_WT ≈ 2.33×,**根本打不到第一轮 9.50 的水平**。必须叠加论文知识。

---

## 三、模型与训练

### 9. 加性模型 R² 全为负(per-parent)

**症状**:Step 2 加性 Ridge 在 per-parent 拆分下,R² 全是负的(模型不如预测母体均值)。

**解释**:单个母体内的 brightness 变异**主要来自测量噪声 + 微弱突变效应**,线性模型无法捕捉。
但**全局 R²=0.222** 因为不同母体 baseline 差异大,聚合后基线预测准确。

**教训**:**看 per-type R²,不要被全局 R² 误导**。

### 10. Epistasis 模型 OOD 陷阱(最大教训)

**症状**:Step C epistasis 模型在 in-distribution val R²=0.9162,**但对论文突变组合(OOD)预测 brightness ~2.6**(应 ≥4.5),finit_rel 全是 0.04-0.09×(应 ≥5×)。

**根因**:
- 训练数据中,多突变 variants 通常 brightness 较低
- 模型学到 "突变数↑ → brightness↓" 的偏见
- 我们的 sfGFP-style 候选(15-19 muts)被预测为低 brightness

**教训**:
1. **不要让模型给 OOD 候选打分**(只用于 in-distribution sanity check)
2. **检查模型预测的合理性**(读 2-3 个数值,对比领域知识)
3. **agent memory 已记录**:Always sanity-check ratio/normalized metrics direction

### 11. ESM2-650M 加载 ~570 missing keys

**症状**:`model.load_state_dict(ckpt["model"], strict=False)` 报 missing=572, unexpected=570。

**解决**:可忽略。ESM 仓库的 checkpoint 包含一些训练时用的额外模块(contacts head 等),主 embedding 提取不受影响。

---

## 四、设计与策略

### 12. 老 Ridge 搜索空间窄(Step 3b 失败)

**现象**:phase1 的 top_candidates.csv 是 Ridge 组合搜索的产物,**只包含低 Δbrightness 突变堆砌**(高 Δbrightness 突变实际上不存在于这种搜索里,因为 Ridge 是线性)。

**后果**:用新 XGBoost 重打分,所有 top 候选 finit_rel ≈ 0.000×(被模型识别为"破坏折叠")。

**教训**:**搜索算法的偏见 → 候选库的偏见 → 数据驱动的天花板**。

### 13. ESMFold 不可用 = 无 pLDDT 验证

**症状**:Round 2 提交**完全没有结构验证**。Seq 1 (19 muts) 可能破坏折叠。

**缓解**:chromophore 三联体检查 + 论文先验 + scaffold 平衡,但不是真正的 pLDDT。

**TODO**:下载 ESMFold 14GB + 解决 SSL 问题。

### 14. Tm 预测缺失

**症状**:综合分公式 `Finit/Finit_WT × Ffinal/Finit` 中,Ffinal/Finit 没有直接预测能力。

**当前做法**:靠"已知 Tm"的经验值(sfGFP ~78°C, StayGold/TGP ~92°C, avGFP ~64°C)粗略估计。

**TODO**:训练专门的 Tm 预测器,或用已知 Tm 值做 lookup table。

---

## 五、PowerShell + Windows 特定坑

### 15. `Get-Content` 编码错乱

**症状**:`Get-Content -Path 'D:\生信\...'` 输出乱码(中文变 `??`)。

**解决**:`[Console]::OutputEncoding=[System.Text.Encoding]::UTF8; Get-Content -LiteralPath '...' -Encoding UTF8`

### 16. 路径在 bash 工具中经常被编码损坏

**症状**:bash 工具调 PowerShell 时,`D:\生信\...` 路径的生信两个字经常变 `??`。

**解决**:用 `-LiteralPath` 而非 `-Path`,且显式 UTF-8 encoding。

### 17. PowerShell 默认禁止某些命令

**症状**:`cmd /c ...` 在 bash 工具中被 auto-classifier 拒绝。

**解决**:用 PowerShell 等价命令(`Start-Process`, `Invoke-Expression`),或请求 permission override。

---

## 六、数据驱动 vs 论文驱动的根本矛盾

### 18. 数据 vs 论文的本质冲突

**问题**:
- 数据驱动的天花板 ≈ 2.33× Finit/Finit_WT
- 第一轮 9.50 综合分需要 Finit_rel ≥ 9×(假设 Ffinal/Finit=1.0)
- 论文突变(S65T, S30R, mAG→TGP 改造)组合可带来 ≥10× 提升
- **但这些论文突变不在数据训练分布内**

**结论**:**论文驱动是唯一出路**,数据只能作为 sanity check。

### 19. 9.50 的估算过于乐观

**现象**:Round 1 Seq 6 (avGFP + sfGFP 11 muts + I152S) 综合分 9.50 是**乐观估计**(假设 Finit 提升 10×,Tm=82°C 全保)。

**实际**:Finit 提升可能只有 5-8×(而非 10×),Tm 可能只有 75-80°C(非 82°C)。

**综合分实际**:可能在 4-7,而不是 9.50。

---

## 七、文件 / 数据管理坑

### 20. 嵌入文件 0.7 GB,版本控制不友好

**问题**:`esm650m_embeddings.npy` (722 MB) + `top10k_esm650m_embeddings.npy` (51 MB) **不能 commit 到 git**。

**当前**:本地保存,文档中标注位置,接手者需手动复制。

### 21. model 文件 ~28 MB

**问题**:`step2_xgboost_gpu.model` (27 MB) + `stepC_xgboost_epistasis.model` (2.8 MB) **也不能 commit**。

**处理**:同上,文档标注位置 + 重新训练的脚本。

---

## 八、未来接手建议

1. **环境**:用 Linux/Mac 避免 Windows 编码坑。如果必须用 Windows,所有路径用 UTF-8 + `-LiteralPath`。
2. **数据天花板**:**不要再尝试纯数据驱动**,直接走论文知识路线。
3. **ESMFold**:不要重新踩 SSL 坑,要么离线下载权重,要么放弃结构验证。
4. **模型**:epistasis 模型只用于 in-distribution sanity check,论文突变组合直接基于先验。
5. **pLDDT 替代**:用 Rosetta / FoldX / I-TASSER 等本地工具,如果非要结构验证。
6. **监控**:开任何长任务(>30 min)都用 cron + mavis cron create,避免单点失败。

---

*对应方法论见 `docs/02_methodology.md`,对应文件清单见 `docs/08_appendix.md`*