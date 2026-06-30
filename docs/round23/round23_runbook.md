# Round 23 Runbook

## 复现步骤

### Step 1: 上传脚本
```bash
gssh cp "D:\workspace\r23_server.py" 9ca7acb1b94c:/root/autodl-tmp/r23_server.py
```

### Step 2: 启动
```bash
gssh run 9ca7acb1b94c "cd /root/autodl-tmp && stdbuf -oL -eL python3 -u r23_server.py > /tmp/r23.log 2>&1"
```

### Step 3: 监控
```bash
gssh --json exec 9ca7acb1b94c "tail -30 /tmp/r23.log"
```

### Step 4: 下载结果
```bash
mkdir -p D:/workspace/round23
gssh cp 9ca7acb1b94c:/root/autodl-tmp/r23/submission_r23.csv D:/workspace/round23/submission_r23.csv
gssh cp 9ca7acb1b94c:/root/autodl-tmp/r23/final_6_r23.json D:/workspace/round23/final_6_r23.json
```

### Step 5: 分析对比
```bash
python D:/workspace/round23/analyze_r23.py
```

## 关键参数

```python
NUM_SEQ_PER_TEMP = 150
TEMPS = [0.1, 0.3, 0.5, 0.7, 1.0]  # 5 档高温度
FIXED = [1, 65, 66, 67, 96, 222]
RECYCLES_SCREEN = 8
BATCH = 25
```

## 流程

```
3 R20 Top 父代 × 5 温度 × 150 候选 = 2250 候选
       ↓
ESMFold r=8 评估
       ↓
1201 passed (53.4%)
       ↓
Sort → Top 6
       ↓
R23 Top 1 = 0.9419
```

## 预计总时间

~5.3 小时 (单任务在 A800 上)

---

*最后更新: 2026-06-29*
