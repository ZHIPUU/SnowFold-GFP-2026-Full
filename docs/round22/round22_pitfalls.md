# Round 22 踩坑记录

## 坑 1: MPNN fixed_positions key (Windows 兼容)

**症状**: 任务 `e1653c347366` exit code 137 (killed)

**根因**: MPNN 在 Windows 下用 `biounit.rfind("/")` 找路径分隔符。Windows 路径 `D:\workspace\r22\pdbs\r22_p1.pdb` 用 `\` 而非 `/`，所以 `rfind("/")` 返回 `-1`，整个路径被当成 PDB name。

**修复**: 把 PDB 路径中的 `\` 替换为 `/`：
```python
pdb_path_mpnn = pdb_path.replace("\\", "/")
```

**教训**: 任何调用 ProteinMPNN 的脚本在 Windows 上**必须**做 path 转换。MPNN 内部不规范化路径。

---

## 坑 2: 输出缓冲（buffering）

**症状**: 日志卡在 `700/1000` 十几分钟不更新,但 GPU 100% 在跑。

**根因**: `gssh run` 通过 `tee` 重定向 stdout 到日志文件，pipe 模式下 Python stdout 变为全缓冲（即使 `flush=True`）。

**修复**: 用 `stdbuf -oL` 强制 line buffering：
```bash
stdbuf -oL -eL python3 -u r22_long.py > /tmp/r22.log 2>&1
```

**教训**: 长任务必须用 `stdbuf -oL` 强制 line-buffered 输出，避免缓冲假象。

---

## 坑 3: gssh cp PowerShell 引号

**症状**: `gssh cp D:\workspace\file.py session:/path/` 报"一端必须是 `<session-id>:<path>`，另一端是本地路径"。

**根因**: gssh cp 对 PowerShell quote 处理敏感。用 `D:/workspace` 而非 `D:\workspace`。

**修复**: 统一用 `D:/workspace/file.py` 风格 + 双引号。

---

## 坑 4: MemoryMPNN subprocess stderr 截断

**症状**: `subprocess.run(capture_output=True)` 拿不到 stderr 详细错误。

**修复**: 检测 `seqs/*.fa` 文件存在性来判断 MPNN 成功与否。

---

## 坑 5: ESM2-650M / ThermoMPNN 下载失败

**症状**: HuggingFace 镜像 (`hf-mirror.com`, `ghproxy.com`, `gh-proxy.polaris-lab.com`) 全部失败。

**根因**: 网络限制。

**影响**: ThermoMPNN ΔΔG、ESM2-650M Fine-tune、ESM3 GFP 生成、RFdiffusion3 全部未跑成。

**教训**: 长跑任务前先用 `curl -I` 测试镜像可用性。

---

## 坑 6: MPNN fixed_positions key 名称匹配

**症状**: `KeyError: 'round23\\pdbs\\r23_p1'`

**根因**: PDB 路径含目录前缀，MPNN 内部用 `biounit.rfind("/")` 提取 name。

**修复**: `pdb_key_simple = os.path.basename(pdb_path).replace(".pdb", "")`，且 PDB 路径必须 forward slash。

---

*记录者: Trae AI Agent (Claude) | 时间: 2026-06-29*
