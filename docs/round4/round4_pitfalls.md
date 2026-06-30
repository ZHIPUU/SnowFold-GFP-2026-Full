# Round 4 踩坑指南 (Pitfalls & Solutions)

> **目的**: 详细记录 Round 4 遇到的所有问题, 帮 Round 5 避免重蹈覆辙
> **每个坑**: 症状 → 根因 → 解决 → 教训

---

## 🕳️ 坑 #1: PyTorch 装成 CPU 版 (最严重)

**症状**:
```python
>>> import torch
>>> torch.__version__
'2.12.1+cpu'  # 错! 应该是 +cu128
>>> torch.cuda.is_available()
False  # 致命
```

**根因**:
- 国内直接 `pip install torch` 会装 PyPI 上的 CPU 版
- 即使指定 `--index-url https://download.pytorch.org/whl/cu128`, 国内网络可能超时
- PyPI 不提供 CUDA wheel, 只有 CPU 版

**解决**:
```powershell
# 用上海交大镜像 (推荐)
pip install torch torchvision --index-url https://mirror.sjtu.edu.cn/pytorch-wheels/cu128 --trusted-host mirror.sjtu.edu.cn

# 验证
python -c "import torch; print(torch.__version__, torch.cuda.is_available())"
# 期望: 2.11.0+cu128 True
```

**阿里云的坑**:
```powershell
# 错误 (阿里云 PyTorch-wheels 目录结构不兼容 PEP 503)
pip install torch torchvision --index-url https://mirrors.aliyun.com/pytorch-wheels/cu128
# 报错: ERROR: Could not find a version that satisfies the requirement torch
```

```powershell
# 正确 (阿里云要用 --find-links)
pip install torch torchvision --find-links https://mirrors.aliyun.com/pytorch-wheels/cu128/ --trusted-host mirrors.aliyun.com
```

**教训**:
- 国内 PyTorch 安装, **必须用支持的镜像源**:
  - 上海交大: 完整支持 `--index-url`
  - 阿里云: 只支持 `--find-links`
  - 清华/中科大: 不支持 PyTorch CUDA wheel
- 不确定时, 先到镜像源 URL 看看有无 `cu128/` 子目录

**Round 3 文档错误**:
> "Python 3.14 不能跑模型"

实际: Python 3.14 完全兼容 cu128, **问题在 pip 装了 CPU 版**!

---

## 🕳️ 坑 #2: ProteinMPNN 中文路径致命

**症状**:
```
FileNotFoundError: 'C:\\Temp\\mpnn_work\\out\\T01//seqs/D:\\生信\\2026Protein Design\\work\\round4\\pdbs\\2B3P.fa'
```

**根因**:
- `protein_mpnn_run.py` 把 PDB 完整路径作为 `name` 嵌入输出路径
- 路径含中文时, 字符串拼接或路径解析出错
- `D:\\生信\\...` 被错误处理为 `D:\\閻㈢喍淇\...` (mojibake)

**解决**:
1. **整个 ProteinMPNN 工作区搬到纯英文路径**:
   ```
   C:\Temp\mpnn_work2\
   ├── ProteinMPNN\   # 从项目目录复制源码
   ├── 2B3P.pdb       # PDB 直接放工作目录根
   ├── parsed.jsonl
   └── out_T01\
   ```

2. **强制切换工作目录**:
   ```python
   import os
   os.chdir("C:\\Temp\\mpnn_work2")  # 切到纯英文路径
   ```

3. **所有 subprocess 调用指定 cwd**:
   ```python
   subprocess.run(cmd, cwd="C:\\Temp\\mpnn_work2")
   ```

**教训**:
- 任何命令行工具 (ProteinMPNN, ESM, PyRosetta 等) 都可能受中文路径影响
- 解决方案: 整个工具链放在 `C:\Temp\xxx\` 纯英文目录
- 项目目录用中文没事, 只在调用外部工具时换路径

**坑中坑**:
- 即使 PDB 在 `D:\生信\2026Protein Design\work\round4\pdbs\2B3P.pdb`, ProteinMPNN 也能跑 (因为我们 cd 到 mpnn_work2 之前已复制了 PDB)
- 关键是 **ProteinMPNN 的输出路径在 cwd 下**, 不能含中文

---

## 🕳️ 坑 #3: ProteinMPNN X 占位符陷阱

**症状**:
```python
# 期望: 实际 MPNN 设计的序列
seq = "SKGEELFTGVVPIK...TYG..."  # 含完整 TYG
# 实际: MPNN 输出
seq = "SKGEELF...TLXXXVQCF..."  # 65-67 是 XXX!
```

**根因**:
- ProteinMPNN 只输出"非固定位"的序列
- 固定位的字符是 `X` (占位符, 表示"这是被锁定的位")
- 我们必须用 PDB 原始序列填回 X

**解决**:
```python
# 读取 PDB 原始序列 (parsed.jsonl)
with open("parsed.jsonl") as f:
    p = json.loads(f.readline())
pdb_wt_with_x = p["seq_chain_A"]  # 包含 X

# 找 PDB 与 sfgfp_wt 的对齐位置
for offset in range(0, 15):
    matches = sum(1 for a, b in zip(pdb_wt_with_x[:30], sfgfp_wt[offset:offset+30]) if a == b or a == 'X')
    if matches >= 28: break

# 用 sfgfp_wt 填回 X
wt_aligned = sfgfp_wt[offset+1:offset+1+len(pdb_wt_with_x)]
for mpnn_seq in design_seqs:
    filled = "".join(wt_aligned[i] if c == 'X' else c for i, c in enumerate(mpnn_seq))
    full_seq = sfgfp_wt[:offset+1] + filled + sfgfp_wt[offset+1+len(filled):]
```

**教训**:
- ProteinMPNN 输出需要后处理 (填回 X + 补 N 端)
- 直接丢弃有 X 的序列 = 100% 错 (大部分是 X)
- **正确处理后才发现 64 条 MPNN 设计序列, 全部保留**

---

## 🕳️ 坑 #4: mBaoJin PDB chromophore 缺失

**症状**:
```python
>>> seq = "ENMASTPFKFQLKGTINGKSFTVEGEGEGNSHEGSHKGKYVCTSGKLPMSWAALGTTF-MKYYTKYPSGL..."
>>> # 期望 GYG 在某处
>>> "GYG" in seq
False  # 找不到!
>>> # 检查 pos 35-40
>>> seq[35:40]
'F-MKY'  # 实际是 'F-MKYY', 不含 GYG
```

**根因**:
- PDB 8QBJ (mBaoJin) 的晶体结构中, chromophore 已自催化成熟
- 序列对应位置的残基被破坏/重排
- 这是 X-ray 解析的"成熟态" 特征, 不是 bug

**解决**:
- mBaoJin 跳过 ProteinMPNN, 用 ESMFold 评估已知 mBaoJin 序列
- mBaoJin WT + 已知突变 (mBaoJin paper 的 8 突变) 即可
- 接受 ESMFold pLDDT 偏低 (StayGold 训练样本少)
- 文献 Tm 92°C 是直接证据, 不依赖 ESMFold

**教训**:
- 某些蛋白的晶体结构"不适合" 直接做反折叠 (成熟 chromophore 破坏)
- 退路: 用 ESMFold 直接评估 + 文献 Tm 作为硬证据

---

## 🕳️ 坑 #5: avGFP 2WUR 起始残基偏移

**症状**:
```python
# 期望: avGFP chromophore TYG 在 pos 65-67 (1-based, 与 sfGFP 同步)
fixed_pos = "65 66 67 96 145 148 167 203 205 222"  # 按文献编号

# 实际: 2WUR chain A 起始 KGE (avGFP WT 第 3 位)
#       chromophore TYG 在 1-based pos 36-38!
```

**根因**:
- 不同 PDB 结构的 residue 编号系统不一致
- 2WUR 起始 K (avGFP pos 3), 但 2WUR 内部 pos 1 已经是 K
- 我们的 fixed_pos 错了, **固定了非 chromophore 位置**

**诊断**:
```python
# 解析后先看 chain A 实际序列
with open("parsed.jsonl") as f:
    p = json.loads(f.readline())
seq = p["seq_chain_A"]
# 找 chromophore
idx = seq.find("TYG")
print(f"TYG at 1-based pos {idx+1}")
# 2WUR: 输出 36 ← 这才是正确 fixed 起点
```

**解决**:
```python
# 用解析后序列定位, 不依赖绝对 pos 编号
chromo_idx = seq.find("TYG")  # 找到 TYG 的 0-based 位置
fixed_pos = [chromo_idx + offset for offset in range(-5, 6)]  # ±5 邻位
# 添加关键功能位 (相对 chromophore 偏移, 不是绝对 pos)
fixed_pos.extend([chromo_idx + 31, chromo_idx + 80, chromo_idx + 138, chromo_idx + 157])
```

**教训**:
- **永远用解析后的真实序列定位功能位**
- **不要套用文献的 pos 编号** (不同 PDB 编号系统差异巨大)
- 关键位可以相对偏移定位 (如 R96 = chromophore + 31)

---

## 🕳️ 坑 #6: ESMFold pLDDT 提取错误

**症状**:
```python
>>> output = model(tokens)
>>> plddt_raw = output.plddt.numpy()[0]  # 错误!
# 实际 plddt shape = (1, L, 37), 不是 (L, 37)
```

**真实形状**:
- `output.plddt` shape = `[1, L, 37]`
- `output.plddt.numpy()[0]` = `(L, 37)`
- 每残基 37 个原子 (atom37 表示)
- 需要 `mean(axis=1)` 得到 per-residue pLDDT

**正确代码**:
```python
plddt_raw = output.plddt.cpu().numpy()[0]  # (L, 37)
atom_mask = output.atom37_atom_exists.cpu().numpy()[0]  # (L, 37)
plddt_scaled = plddt_raw * 100.0
masked_sums = (plddt_scaled * atom_mask).sum(axis=1)
masked_counts = atom_mask.sum(axis=1).astype(float)
masked_counts[masked_counts == 0] = 1
plddt_per_res = masked_sums / masked_counts
mean_plddt = float(plddt_per_res.mean())
```

**坑 6.2**: 没用 atom_mask
```python
# ❌ 错
plddt_per_res = plddt_scaled.mean(axis=1)  # 没考虑 N 端缺原子

# ✅ 对
plddt_per_res = (plddt_scaled * atom_mask).sum(axis=1) / atom_mask.sum(axis=1)
```

**坑 6.3**: 标量 ptm 误用 list 索引
```python
# ❌ 错
ptm = output.ptm[0].cpu().numpy()  # IndexError: invalid index of 0-dim tensor

# ✅ 对
ptm = output.ptm.cpu().item()
```

**教训**:
- 总是 `print(output)` 或 `print(type(output.plddt), output.plddt.shape)` 看实际形状
- 不用假设文档说的形状, 实际跑一次确认

---

## 🕳️ 坑 #7: FP16 vs FP32 精度

**症状**:
- 第一次跑 ESMFold, pLDDT 全部 0.3-0.5
- 看起来像"全错"

**根因**:
```python
model.esm = model.esm.half()  # FP16
# FP16 严重降低 pLDDT 精度
```

**解决**:
```python
# 不转换, 保持默认 FP32
model = model.cuda()
# 或显式恢复
model.esm = model.esm.float()
```

**教训**:
- pLDDT 头对 FP16 特别敏感
- 即使显存紧张, 也不要对 ESMFold 用 FP16
- ESM-IF1 同样要 FP32

---

## 🕳️ 坑 #8: PowerShell 路径与中文引号

**症状**:
```powershell
python "work\round4\my_script.py" --name "我的蛋白"
# 可能报错: 系统找不到路径 "D:\生信\..." 或参数被截断
```

**根因**:
- PowerShell 默认编码 GBK, 与 Python UTF-8 不兼容
- 引号嵌套问题 (`"..."` 在 `"..."` 中)
- 多重转义

**解决**:
```powershell
# 永远用脚本文件, 不用 -c 一行命令
python "work\round4\my_script.py"  # 脚本中用 pathlib.Path(r"...") 显式 UTF-8

# 复杂命令
python "D:\生信\2026Protein Design\work\round4\long_script.py"
# ✅ 用绝对路径 + 引号 + .py 扩展名
```

**坑 8.2**: `Find` 命令找不到路径
```powershell
# ❌ 错
Get-ChildItem work/round4
# 可能错 (相对路径 + 引号)

# ✅ 对
Get-ChildItem "D:\生信\2026Protein Design\work\round4"
# 或
Set-Location "D:\生信\2026Protein Design"
Get-ChildItem work/round4
```

**教训**:
- **永远用脚本文件**, 不用 `-c` 一行命令
- Python 脚本用 `pathlib.Path(r"...")` 显式 UTF-8
- PowerShell 路径用**绝对路径**最稳定

---

## 🕳️ 坑 #9: cgreGFP 全面低 pLDDT

**症状**:
```
Y2_cgreGFP_S65T_K163A: pLDDT=30.9
Y3_cgreGFP_R1_combo: pLDDT=30.2
Y1_cgreGFP_S65T: ❌ 命中排除列表
```

**根因**:
- cgreGFP 在 ESMFold 训练集中样本可能较少
- cgreGFP 序列与 avGFP/sfGFP 差异较大 (同源性 ~30%)
- ESMFold 对 cgreGFP 评估置信度低 (pLDDT 30)

**解决**:
- cgreGFP 候选全部剔除
- 用其他骨架 (avGFP, amacGFP, mBaoJin) 代替
- cgreGFP 即使 baseline brightness 高 (4.50), 也比不过高 pLDDT 候选

**教训**:
- pLDDT < 35 = 实际折叠概率 < 60%
- 宁可放弃"理论上亮度高" 的候选, 也要保证"实际能折叠"
- Round 3 错误地保留了 cgreGFP, Round 4 修正

---

## 🕳️ 坑 #10: Tm 估值混淆 (v3 → v4 → v5)

**症状**:
- v3: 所有 MPNN 候选 Tm 估 80°C
- v4: 同 v3
- v5: 按 pLDDT 分级估 72-88°C

**为什么 v3/v4 错?**
- 高 pLDDT (= 68) 反映序列与结构高度自洽
- 物理规律: 高 pLDDT ≈ 天然蛋白 ≈ 高 Tm
- v3/v4 把 MPNN_T01_014 (pLDDT 68) Tm 估 80°C 严重低估

**解决 (v5)**:
```python
if plddt >= 65: tm = 88
elif plddt >= 60: tm = 84
elif plddt >= 55: tm = 81
elif plddt >= 50: tm = 78
else: tm = 72  # 高风险
```

**教训**:
- **多步推理链路中, 任何一步估值偏差都会传导**
- 估值要"上下界合理" (高 pLDDT 不能估 80°C, 至少 85+)
- 写文档时显式说明估值方法, 方便后任校验

---

## 🕳️ 坑 #11: git clone 失败

**症状**:
```
fatal: unable to access 'https://github.com/dauparas/ProteinMPNN.git/': 
Failed to connect to github.com port 443 after 21092 ms: Could not connect
```

**根因**:
- 国内访问 github.com 经常超时/被墙
- 即使能访问, 速度极慢

**解决**:
```bash
# GitHub 镜像 (任选一个)
git clone https://gitcode.com/gh_mirrors/pr/ProteinMPNN.git
# 或
git clone https://gitee.com/mirrors/ProteinMPNN.git
```

**教训**:
- 国内 git clone **必须用镜像**
- gitcode.com 通常最新
- gitee.com 镜像质量参差不齐, 可能找不到

---

## 🕳️ 坑 #12: GPU Out of Memory (OOM)

**症状**:
```
RuntimeError: CUDA out of memory. Tried to allocate 2.00 GiB
```

**根因**:
- ESMFold 推理每序列 ~5GB 显存
- RTX 5080 16GB 最多同时 2-3 序列
- 长序列 (> 400 aa) 更耗

**解决**:
```python
# 1. 减少 chunk size
model.trunk.set_chunk_size(64)  # 默认 128 → 64, 节省显存

# 2. 不在 GPU 上同时保存多个序列
for i, c in enumerate(candidates):
    with torch.no_grad():
        tokens = tokenizer([seq], return_tensors="pt")["input_ids"].cuda()
        output = model(tokens)
    del tokens, output
    torch.cuda.empty_cache()  # 强制释放

# 3. 减少 batch size (目前用 batch=1, 没问题)
```

**教训**:
- 总是 `torch.cuda.empty_cache()` 释放显存
- 16GB GPU 跑 ESMFold 没问题, ProteinMPNN 也没问题
- 多模型同时在 GPU 容易 OOM

---

## 🕳️ 坑 #13: docx/PDF 生成路径问题

**症状**:
- 早期 Round 1-3 用 `md_to_pdf.py` 生成设计文档
- 在 PowerShell 下, 输出文件名含中文可能乱码

**解决**:
- 用 reportlab (Python) 生成 PDF, 不用命令行 pandoc/wkhtmltopdf
- 文件名用英文, 内容用中文
- 输出目录纯英文

**教训**:
- 中文 PDF 用 reportlab, 不用 wkhtmltopdf
- 测试 PDF 在不同查看器中是否正确显示中文

---

## 🕳️ 坑 #14: 排除列表内存爆炸

**症状**:
```python
excl = pd.read_csv("Exclusion_List.csv")  # 135,414 行
excl_seqs = set(excl["Sequence"])
# 内存占用 ~200MB
```

**根因**:
- 13.5 万行 CSV 读入 pandas, 每个 Sequence 字符串 ~250 char
- 字符串对象开销大

**解决**:
```python
# 用 set 查询更快 + 内存优化
excl_seqs = set(pd.read_csv("Exclusion_List.csv")["Sequence"].astype(str).str.strip())
# 查询 O(1), 200MB 内存可接受
```

**教训**:
- 排除列表用 set, 不要 list 查询
- `.str.strip()` 防止末尾空白字符导致漏匹配

---

## 🕳️ 坑 #15: Tm 修正后的乐观场景反而下降

**症状**:
- v3 → v5 修正 Tm 估值后, 乐观场景预测从 2.03 → 1.90 (下降)

**根因**:
- v3 选 Top-6 全部手工设计 (突变少, Tm 高估)
- v5 加入 MPNN (突变多, Tm 保守估)
- 乐观场景下, 手工设计的 Tm 上限高, MPNN 的上限低
- **比较基础不同, 不能直接对比**

**教训**:
- 版本对比时, 还要看"实际比赛风险", 不只是预测分
- v5 在"低风险"上更强 (MPNN 68 pLDDT = 几乎确定折叠)
- v3 在"理论高分"上更强 (全部手工 + Tm 估值乐观)

---

## 🧭 综合建议 (避免 Round 5 重蹈覆辙)

### 优先级排序的"该做"和"不该做"

**✅ 该做**:
1. **提交前 30 分钟**先做一遍完整流程 (避免截止前才发现问题)
2. **所有新工具** 先在小数据集 (1-2 条) 跑通, 再大规模
3. **数据验证**: 读 2-3 行确认格式正确 (尤其 pLDDT、序列)
4. **错误立刻记录**: 别用脑子记, 写文件
5. **保存中间产物**: JSON/CSV, 方便重跑

**❌ 不该做**:
1. **不要信任文档**: 实测才知道 (Round 3 文档错误 3 处)
2. **不要在中文路径跑 MPNN/ESMFold/ESM-IF**
3. **不要给 MPNN 装 CPU 版 PyTorch** (装错用错)
4. **不要固定过多位** (MPNN 设计自由度受限)
5. **不要使用 ML 打分 OOD** (Round 2 教训)
6. **不要盲目追求 v5/v4 之间的预测分** (风险调整后 v5 更优)

### 调试模板

```python
# 标准调试代码
import torch
import json
from pathlib import Path

# 1. 验证环境
print("torch:", torch.__version__)
print("cuda:", torch.cuda.is_available())
if torch.cuda.is_available():
    print("device:", torch.cuda.get_device_name(0))
    print("sm:", torch.cuda.get_device_capability(0))

# 2. 验证数据
ROOT = Path(r"D:\生信\2026Protein Design")
for p in [ROOT / "AAseqs of 5 GFP proteins_20260511.txt",
         ROOT / "Exclusion_List.csv",
         ROOT / "work/round4/submission_round4_v5.csv"]:
    print(f"{p.name}: exists={p.exists()}, size={p.stat().st_size}")

# 3. 验证模型
from transformers import AutoTokenizer, EsmForProteinFolding
try:
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
    print("tokenizer OK")
except Exception as e:
    print(f"tokenizer FAIL: {e}")
```

### 关键文件 (Round 4 维护)

1. **work/round4/submission_round4_v5.csv** - 最终提交, **不要改**
2. **work/round4/final_6_round4_v5.json** - 详细评分, 可改但要重新生成
3. **docs/round4/round4_*.md** - 6 份文档, 接手时必读
4. **work/round4/ProteinMPNN/** - MPNN 源码, 跑了 5 次, 验证可用

### 关键时间点

- 比赛 CFPS 协议 → **联系 root 询问** (5 分钟)
- 设计 PDF → 1-2 小时
- GitHub 仓库 → 30 分钟
- ThermoMPNN 部署 → 30-60 分钟
- PROSS 在线服务 → 30 分钟
- 新跑 MPNN (T=0.5, T=0.7) → 30-60 分钟

**总时间预算**: 假设剩 1 天, 优先做"比赛必交材料", 然后再优化。

---

**完成日期**: 2026-06-22
**作者**: Trae AI Agent
**目的**: 让 Round 5 接手者少走 5 小时的坑
