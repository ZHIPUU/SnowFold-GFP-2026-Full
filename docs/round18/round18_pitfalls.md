# Round 18 踩坑记录

## 坑 1: MPNN fixed_positions 找不到 PDB

**症状**: MPNN 跑出 0 条候选,日志无明显错误

**根因**:
```python
# 错误: pdb_key = "/root/autodl-tmp/r18/pdbs/R14_A_T01_013" (全路径)
# MPNN 内部用 os.path.basename 比较,所以找不到
```

**修复**:
```python
pdb_key = os.path.basename(pdb_path).replace(".pdb", "")
```

## 坑 2: MPNN subprocess stderr 截断

**症状**: `subprocess.run(..., capture_output=True)` 拿不到 stderr

**根因**: MPNN 输出到 stderr 太长,被 subprocess 截断

**解决**: 检测 `glob.glob("seqs/*.fa")` 是否为空,为空则重跑不 capture

## 坑 3: HF ESMFold 必须用 `out.plddt[..., 1]`

**教训**: 不要用 `out.plddt[..., :].mean(-1)`(那是平均 bin 索引,不是真实 pLDDT)

---

*记录者: Trae AI Agent (Claude) | 最后更新: 2026-06-28*
