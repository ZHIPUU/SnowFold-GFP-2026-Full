# Round 19 Runbook

## 复现步骤

### Step 1: 上传并运行
```bash
gssh cp "D:\workspace\r19_pipeline.py" 9ca7acb1b94c:/root/autodl-tmp/r19_pipeline.py
gssh exec 9ca7acb1b94c "cd /root/autodl-tmp && python3 r19_pipeline.py"
```

### Step 2: 下载结果
```bash
gssh cp 9ca7acb1b94c:/root/autodl-tmp/r19/submission_r19.csv "D:\workspace\round19\"
gssh cp 9ca7acb1b94c:/root/autodl-tmp/r19/final_6_r19.json "D:\workspace\round19\"
```

### Step 3: 合规检查
```bash
python "D:\workspace\round19\check_compliance.py"
```

## 关键参数
- 父代: 9 (6 R18 + 3 WT)
- 温度: 0.05 / 0.1 / 0.2 / 0.3 / 0.5
- 候选/温度: 150
- batch_size: 20

## 耗时: ~167 分钟

---

*最后更新: 2026-06-28*
