# Round 16 Runbook

## 复现步骤

### Step 1: 验证 fair-esm vs HuggingFace pLDDT 差异
```python
# 本地 RTX 5080
python work/round16/test_fair_esmfold.py
```

### Step 2: 跑 R17 校准管线
```bash
gssh exec 9ca7acb1b94c "cd /root/autodl-tmp && python3 r17_pipeline.py"
```

### Step 3: 对比输出
- `D:\workspace\round17\submission_r17.csv` — HF 0-1 尺度评分
- R15 旧 `submission_top6.csv` — fair-esm 0-100 尺度(不可对比)

## 关键文件
- [round16_report.md](round16_report.md) — 实验报告
- [round16_pitfalls.md](round16_pitfalls.md) — 踩坑记录
