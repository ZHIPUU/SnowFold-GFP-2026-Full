# Round 20 踩坑记录

## 坑 1: Fixed Position 1 (M) 缺失 ⚠️ **关键**

**症状**: 
- R18/R19 MPNN 输出中 **97% 候选不以 M 开头** (实测 99/250 实测)
- 实际只有 28/1000 候选通过 `seq.startswith("M")` 过滤

**根因**:
```python
# R18/R19 错误写法
FIXED = [65, 66, 67, 96, 222]  # 只固定 chromophore 残基
```

MPNN 默认不固定 position 1，必须显式告诉它 position 1 (M) 不能改变。

**修复**:
```python
FIXED = [1, 65, 66, 67, 96, 222]  # +1 (M) + 5 chromophore
```

**验证**: 修复后 100% 候选以 M 开头，通过率从 ~8% 跃升到 **62.8%**。

**教训**: 任何 GFP 设计必须确保 position 1 (M) 被 MPNN 固定，这是项目规则第一节"首位必须为 M"的硬要求。

---

## 坑 2: print 输出缓冲

**症状**: 即使 `python3 -u` (unbuffered) + `flush=True`，gssh `run` 任务的 stdout 仍可能缓冲到几十分钟才 flush。

**根因**: gssh `run` 用 `tee` 重定向 stdout 到日志文件，pipe 模式下 Python stdout 会变成**全缓冲**（不是行缓冲）。

**修复**:
```bash
stdbuf -oL -eL python3 -u r22_long.py > /tmp/r22.log 2>&1
```

`stdbuf -oL` 强制 line-buffered 输出。配合 `flush=True` 双重保险。

**教训**: 长任务监控时不能只看 stdout，要同时查：
- `ps aux | grep python`
- `nvidia-smi --query-gpu=utilization.gpu`
- `/proc/<PID>/io` (write_bytes)

---

## 坑 3: gssh cp Windows 路径

**症状**: 
```bash
gssh cp "D:\workspace\r17_pipeline.py" 9ca7acb1b94c:/root/autodl-tmp/
```
错误信息: "一端必须是 `<session-id>:<path>`，另一端是本地路径"

**修复**:
- 必须**双引号**包路径
- 用 `D:/workspace/` 风格，不用 `D:\workspace\` 反斜杠

**教训**: gssh cp 对 PowerShell quote 处理敏感，统一用 `D:/...` 风格。

---

## 坑 4: ProteinMPNN subprocess stderr 截断

**症状**: `subprocess.run(capture_output=True)` 拿不到 stderr，或者 stderr 被截断。

**修复**:
```python
files = list_fa_files(outdir)  # 检查输出文件
if not files:
    subprocess.run(cmd, timeout=600)  # 重跑不 capture
```

**教训**: 用**输出文件存在性**判断 MPNN 成功与否，比 stderr 更可靠。

---

## 坑 5: MPNN fixed_positions key

**症状**: MPNN 输出 1 个空 .fa 文件 (0 候选)。

**修复**:
```python
# 错误
pdb_key = pdb_path.replace(".pdb", "")  # 全路径, MPNN 找不到

# 正确
pdb_key = os.path.basename(pdb_path).replace(".pdb", "")
```

**教训**: `--fixed_positions_jsonl` 的 key 必须匹配 MPNN internal PDB basename。

---

## 坑 6: glob() 解析名字带逗号

**症状**: `glob.glob("seqs/*.fa")` 返回空，但 `os.walk(seqs_dir)` 能找到。

**根因**: MPNN 输出目录名包含 R19 父代的 `T=0.2, sample=25, ...` 字符串，有逗号。glob 在某些 shell 下解析失败。

**修复**: 用 `os.walk()` 替代 `glob.glob()`。

**教训**: 名字带特殊字符 (`,`, ` `, `&`) 用 `os.walk` 永远安全。

---

*记录者: Trae AI Agent (Claude) | 最后更新: 2026-06-29*
