# Round 17 踩坑记录

## 坑 1: openfold 编译失败 (DLL 错误)

**症状**: `pip install openfold` 报错 `git returned 4294967290`

**根因**: openfold 需要从源码编译 C++/CUDA 扩展,Windows MSVC + git 路径冲突

**教训**: **完全不需要 openfold**, HuggingFace ESMFold 已经是生产级实现

## 坑 2: git clone github 失败

**症状**: `git clone https://github.com/...` 超时

**根因**: A800 平台未启用 network turbo

**解决**: `source /etc/network_turbo`,或换用 ghproxy.com 镜像

## 坑 3: ProteinMPNN pip install 失败

**症状**: `pip install git+https://github.com/dauparas/ProteinMPNN.git` 超时

**解决**: **直接本地 cp 核心文件**:
- `protein_mpnn_run.py`
- `protein_mpnn_utils.py`
- `vanilla_model_weights/v_48_020.pt`

## 坑 4: MPNN fixed_positions key 错误

**症状**: MPNN 跑出 0 条候选

**根因**: `--fixed_positions_jsonl` 的 key 必须是 PDB **basename**,不是全路径

**修复**:
```python
pdb_key = os.path.basename(pdb_path).replace(".pdb", "")
# 而不是
pdb_key = pdb_path.replace(".pdb", "")
```

---

*记录者: Trae AI Agent (Claude) | 最后更新: 2026-06-28*
