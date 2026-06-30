# Round 17 Runbook

## 复现步骤

### Step 1: 上传脚本
```bash
gssh cp "D:\workspace\r17_pipeline.py" 9ca7acb1b94c:/root/autodl-tmp/r17_pipeline.py
```

### Step 2: 跑 R17
```bash
gssh exec 9ca7acb1b94c "cd /root/autodl-tmp && python3 r17_pipeline.py"
```

### Step 3: 下载结果
```bash
gssh cp 9ca7acb1b94c:/root/autodl-tmp/r17/submission_r17.csv "D:\workspace\round17\"
```

## 关键修复
- `predict()` 函数用 `out.plddt[..., 1]` 直接 0-1 尺度
- 评分公式直接 `0.3*pLDDT`,不再 `/100`

## 已知 bug 修复
- R15 fair-esm `out["plddt"]` 输出 0-100 整数
- R17 HF `out.plddt[..., 1]` 输出 0-1 浮点
- **两者相差 100 倍**,不能混用

---

*最后更新: 2026-06-28*
